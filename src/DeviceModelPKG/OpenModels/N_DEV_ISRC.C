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
// Purpose        : Independent current source
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
#include <N_DEV_ISRC.h>
#include <N_DEV_MatrixLoadData.h>
#include <N_DEV_SolverState.h>
#include <N_DEV_SourceData.h>
#include <N_DEV_Message.h>
#include <N_ERH_ErrorMgr.h>

#include <N_LAS_Vector.h>
#include <N_LAS_Matrix.h>
#include <N_UTL_FeatureTest.h>
#include <N_UTL_MachDepParams.h>


#include <N_UTL_Math.h>

namespace Xyce {
namespace Device {
namespace ISRC {

void Traits::loadInstanceParameters(ParametricData<ISRC::Instance> &p)
{
  // DC parameters
  p.addPar ("DCV0", 0.0, &ISRC::Instance::DCV0)
    .setOriginalValueStored(true)
    .setUnit(U_VOLT)
    .setDescription("DC Current");

  // Pulse parameters
  p.addPar ("V0", 0.0, &ISRC::Instance::par0)
    .setUnit(U_AMP)
    .setDescription("Offset Current");

  p.addPar ("V1", 0.0, &ISRC::Instance::par0)
    .setUnit(U_AMP)
    .setDescription("Initial Current");

  p.addPar ("V2", 0.0, &ISRC::Instance::par1)
    .setUnit(U_AMP)
    .setDescription("Pulsed Current");

  p.addPar ("TD", 0.0, &ISRC::Instance::par2)
    .setUnit(U_SECOND)
    .setDescription("Delay");

  p.addPar ("TR", 0.0, &ISRC::Instance::par3)
    .setUnit(U_SECOND)
    .setDescription("Rise Time");

  p.addPar ("TF", 0.0, &ISRC::Instance::par4)
    .setUnit(U_SECOND)
    .setDescription("Fall Time");

  p.addPar ("PW", 0.0, &ISRC::Instance::par5)
    .setUnit(U_SECOND)
    .setDescription("Pulse Width");

  p.addPar ("PER", 0.0, &ISRC::Instance::par6)
    .setUnit(U_SECOND)
    .setDescription("Period");

  // Sin parameters
  p.addPar ("VA", 0.0, &ISRC::Instance::par1)
    .setUnit(U_AMP)
    .setDescription("Amplitude");

  p.addPar ("FREQ", 0.0, &ISRC::Instance::par3)
    .setUnit(U_SECM1)
    .setDescription("Frequency");

  p.addPar ("THETA", 0.0, &ISRC::Instance::par4)
    .setDescription("Theta");

  p.addPar ("PHASE", 0.0, &ISRC::Instance::par5)
    .setDescription("Phase");

  // Exp parameters
  p.addPar ("TD1", 0.0, &ISRC::Instance::par2)
    .setUnit(U_SECOND)
    .setDescription("Rise Delay Time");

  p.addPar ("TAU1", 0.0, &ISRC::Instance::par3)
    .setUnit(U_SECOND)
    .setDescription("Rise Time Constant");

  p.addPar ("TD2", 0.0, &ISRC::Instance::par4)
    .setUnit(U_SECOND)
    .setDescription("Fall Delay Time");

  p.addPar ("TAU2", 0.0, &ISRC::Instance::par5)
    .setUnit(U_SECOND)
    .setDescription("Fall Time Constant");

  // AC parameters
  p.addPar ("ACMAG", 0.0, &ISRC::Instance::ACMAG)
    .setUnit(U_VOLT)
    .setDescription("Amplitude");

  p.addPar ("ACPHASE", 0.0, &ISRC::Instance::ACPHASE)
    .setDescription("Phase");

  // SFFM parameters
  p.addPar ("FC", 0.0, &ISRC::Instance::par2)
    .setUnit(U_SECM1)
    .setDescription("Carrier Frequency");

  p.addPar ("FS", 0.0, &ISRC::Instance::par4)
    .setUnit(U_SECM1)
    .setDescription("Signal Frequency");

  p.addPar ("MDI", 0.0, &ISRC::Instance::par3)
    .setDescription("Modulation Index");

  // PWL params
  p.addPar ("R", 0.0, &ISRC::Instance::REPEATTIME)
    .setUnit(U_SECOND)
    .setDescription("Repeat Time");

  p.addPar ("T", 0.0, &ISRC::Instance::T)
    .setUnit(U_SECOND)
    .setDescription("Time");  // time-voltage pairs

  p.addPar ("V", 0.0, &ISRC::Instance::V)
    .setUnit(U_AMP)
    .setDescription("Current"); // time-voltage pairs

  // PAT params
  p.addPar ("VHI", 0.0, &ISRC::Instance::par0)
   .setUnit(U_VOLT)
   .setCategory(CAT_NONE)
   .setDescription("High Voltage Value");

  p.addPar ("VLO", 0.0, &ISRC::Instance::par1)
   .setUnit(U_VOLT)
   .setCategory(CAT_NONE)
   .setDescription("Low Voltage Value");

  p.addPar ("TSAMPLE", 0.0, &ISRC::Instance::par5)
   .setUnit(U_SECOND)
   .setCategory(CAT_NONE)
   .setDescription("Pulse Width");

  p.addPar ("DATA", "", &ISRC::Instance::DATA)
   .setUnit(U_NONE)
   .setCategory(CAT_NONE)
   .setDescription("Data Pattern");

  p.addPar ("RB", 1, &ISRC::Instance::RB)
   .setUnit(U_NONE)
   .setCategory(CAT_NONE)
   .setDescription("Starting Bit When Repeating");

  // Set up non-double precision variables:
  p.addPar ("TRANSIENTSOURCETYPE", (int)_DC_DATA, &ISRC::Instance::TRANSIENTSOURCETYPE)
    .setGivenMember(&ISRC::Instance::TRANSIENTSOURCETYPEgiven);

  p.addPar ("ACSOURCETYPE", (int) _AC_DATA, &ISRC::Instance::ACSOURCETYPE)
    .setGivenMember(&ISRC::Instance::ACSOURCETYPEgiven);

  p.addPar ("DCSOURCETYPE", (int) _DC_DATA, &ISRC::Instance::DCSOURCETYPE)
    .setGivenMember(&ISRC::Instance::DCSOURCETYPEgiven);

  p.addPar ("NUM", 0, &ISRC::Instance::NUM);
}

void Traits::loadModelParameters(ParametricData<ISRC::Model> &p)
{
}


std::vector< std::vector<int> > Instance::jacStamp;
std::vector< std::vector<int> > Instance::jacStampPDE;

//-----------------------------------------------------------------------------
// Function      : Instance::Instance
// Purpose       : instance block constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/06/00
//-----------------------------------------------------------------------------
Instance::Instance(
  const Configuration & configuration,
  const InstanceBlock & instance_block,
  Model &               model,
  const FactoryBlock &  factory_block)
  : SourceInstance(instance_block, configuration.getInstanceParameters(), factory_block),
    model_(model),
    li_Pos(-1),
    li_Neg(-1),
    li_branch_data(0),
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
    ACMAG(1.0),
    ACPHASE(0.0),
    firstTimeload(true),
    freqVarsLoaded(false)
{
  numIntVars = 0;
  numExtVars = 2;
  numStateVars = 0;
  setNumBranchDataVars(0);             // by default don't allocate space in branch vectors
  numBranchDataVarsIfAllocated = 1;    // this is the space to allocate if lead current or power is needed.

  if( jacStamp.empty() )
  {
    jacStamp.resize(2);

    // PDE supporting stamp.  This includes diagonal elements, needed by the
    // 2-level Newton.
    jacStampPDE.resize(2);
    jacStampPDE[0].resize(1);
    jacStampPDE[1].resize(1);
    jacStampPDE[0][0] = 0;
    jacStampPDE[1][0] = 1;
  }

  // Set params to constant default values:
  setDefaultParams ();

  // Set params according to instance line and constant defaults from metadata:
  setParams (instance_block.params);

  const SolverState &solver_state = factory_block.solverState_;
  const DeviceOptions &device_options = factory_block.deviceOptions_;

  // Set any non-constant parameter defaults:
  if (ACSpecified_ && ACSOURCETYPEgiven)
  {
    acSourceData_ = new ACData(*this, instance_block.params, solver_state, device_options);
  }


  dcSourceData_ = new ConstData(*this, instance_block.params, solver_state, device_options);

  if (HBSpecified_ || TRANSIENTSOURCETYPEgiven)
  {
    switch (TRANSIENTSOURCETYPE)
    {
      case _SIN_DATA:
        tranSourceData_ = new SinData(*this, instance_block.params, solver_state, device_options);
        break;

      case _EXP_DATA:
        tranSourceData_ = new ExpData(*this, instance_block.params, solver_state, device_options);
        break;

      case _PULSE_DATA:
        tranSourceData_ = new PulseData(*this, instance_block.params, solver_state, device_options);
        break;

      case _PWL_DATA:
        tranSourceData_ = new PWLinData(*this, instance_block.params, solver_state, device_options);
        break;

      case _PAT_DATA:
        tranSourceData_ = new PatData(*this, instance_block.params, solver_state, device_options);
        break;

      case _SFFM_DATA:
        tranSourceData_ = new SFFMData(*this, instance_block.params, solver_state, device_options);
        break;

      case _DC_DATA:
        tranSourceData_ = 0; // this forces us to use the DC source pointer
        break;

      default:
        UserError(*this) << "Cannot identify source data type for " << getName();
        break;
    }
  }

  processParams();

  // Calculate any parameters specified as expressions:

  updateDependentParameters();
  processParams();
}

//-----------------------------------------------------------------------------
// Function      : Instance::processParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/06/00
//-----------------------------------------------------------------------------
bool Instance::processParams ()
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

// additional Declarations
//-----------------------------------------------------------------------------
// Function      : Instance::registerLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 6/20/02
//-----------------------------------------------------------------------------
void Instance::registerLIDs( const std::vector<int> & intLIDVecRef,
	                               const std::vector<int> & extLIDVecRef)
{
  AssertLIDs(intLIDVecRef.size() == numIntVars);
  AssertLIDs(extLIDVecRef.size() == numExtVars);

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) )
  {
    Xyce::dout() << std::endl << section_divider << std::endl;
    Xyce::dout() << "  ISRCInstance::registerLIDs" << std::endl;
    Xyce::dout() << "  name = " << getName() << std::endl;
  }

  intLIDVec = intLIDVecRef;
  extLIDVec = extLIDVecRef;

  // Now use these lists to obtain the indices into the
  // linear algebra entities.  This assumes an order.    Note that
  // for a current  source, there  will be no Jacobian entries.

  li_Pos = extLIDVec[0];
  li_Neg = extLIDVec[1];

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) )
  {
    Xyce::dout() << "  li_Pos = " << li_Pos << std::endl;
    Xyce::dout() << "  li_Neg = " << li_Neg << std::endl;
    Xyce::dout() << section_divider << std::endl;
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
void Instance::registerStateLIDs(const std::vector<int> & staLIDVecRef )
{
  AssertLIDs(staLIDVecRef.size() == numStateVars);
}

//-----------------------------------------------------------------------------
// Function      : Instance::registerStoreLIDs
// Purpose       : One store var for device current.
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical Systems Modeling
// Creation Date : 03/27/2013
//-----------------------------------------------------------------------------
void Instance::registerStoreLIDs(const std::vector<int> & stoLIDVecRef )
{
}

//-----------------------------------------------------------------------------
// Function      : Instance::registerBranchDataLIDs
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
  if (loadLeadCurrent)
  {
    addBranchDataNode( symbol_table, li_branch_data, getName(), "BRANCH_D");
  }
}


//-----------------------------------------------------------------------------
// Function      : Instance::jacobianStamp
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 9/5/2
//-----------------------------------------------------------------------------
const std::vector< std::vector<int> > & Instance::jacobianStamp() const
{
  if (getSolverState().isPDESystem_)
    return jacStampPDE;

  return jacStamp;
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

	if ( ( fabs(frequency - freq * fIdx) < (frequency * tol  + tol ) ) &&  ( 2* fIdx + 1 <= size_ ) )
	  tmpVal = std::complex<double> ( ftOutData_[ 2* fIdx]/size_ , ftOutData_[ 2* fIdx + 1 ]/size_);

      }
      break;

      case _DC_DATA:
      {
	if (frequency == 0.0 )
	  tmpVal = std::complex<double> ( v0, 0 );

      }
      break;

      default:
        UserFatal(*this) << "Cannot identify source data type for " << getName();
        break;
    }
     // Add RHS vector element for the positive circuit node KCL equ.
    tmpEntry.val = -tmpVal;
    tmpEntry.lid = li_Pos;
    bVec.push_back(tmpEntry);

    tmpEntry.val = tmpVal;
    tmpEntry.lid = li_Neg;
    bVec.push_back(tmpEntry);
  }        
     
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::calculateFDVars
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
// Function      : Instance::loadBVectorsforAC
//
// Purpose       : Loads the B-vector contributions for a single
//                 isrc instance.
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
    acSourceData_->setRealFlag(true);

    acSourceData_->updateSource();
    double source = acSourceData_->returnSource();

    bVecReal[li_Pos] -= source;
    bVecReal[li_Neg] += source;

    acSourceData_->setRealFlag(false);

    acSourceData_->updateSource();
    source = acSourceData_->returnSource();

    bVecImag[li_Pos] -= source;
    bVecImag[li_Neg] += source;
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEBVector
// Purpose       : Loads the F-vector contributions for a single
//                 ISRC instance.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/05/04
//-----------------------------------------------------------------------------
bool Instance::loadDAEBVector ()
{
  double * bVec = extData.daeBVectorRawPtr;

  // get the source value:
  SourceData *dataPtr = dcSourceData_; // by default assume the DC value.
  if ((HBSpecified_ || ( getSolverState().tranopFlag && (!getSolverState().locaEnabledFlag || !DCSOURCETYPEgiven ) ) || getSolverState().transientFlag || (ACSpecified_ && !DCSOURCETYPEgiven ) ) && tranSourceData_ != 0 )
  {
    dataPtr = tranSourceData_;
  }

  double source = 0.0;
  if (dataPtr != 0)
  {
    source = dataPtr->returnSource();
  }
  bVec[li_Pos] -= source;
  bVec[li_Neg] += source;

  if( loadLeadCurrent )
  {
    double * leadF = extData.nextLeadCurrFCompRawPtr;
    double * junctionV = extData.nextJunctionVCompRawPtr;
    double * solVec = extData.nextSolVectorRawPtr;
    leadF[li_branch_data] = source;
    junctionV[li_branch_data] = solVec[li_Pos] - solVec[li_Neg];
  }

  return true;
}

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
  : DeviceModel(MB, configuration.getModelParameters(), factory_block)
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

// additional Declarations

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
  if (DEBUG_DEVICE)
  {
    std::vector<Instance*>::const_iterator iter;
  std::vector<Instance*>::const_iterator first = instanceContainer.begin();
  std::vector<Instance*>::const_iterator last  = instanceContainer.end();

  int i;
  os << std::endl;
  os << "    name     modelName  Parameters" << std::endl;
  for (i=0, iter=first; iter!=last; ++iter, ++i)
  {
    os << "  " << i << ": " << (*iter)->getName() << "      ";
    os << getName();
    os << std::endl;
    if ( (*iter)->tranSourceData_ != 0 )
    {
      (*iter)->tranSourceData_->printOutParams ();
    }

    if ( (*iter)->dcSourceData_ != 0 )
    {
      (*iter)->dcSourceData_->printOutParams ();
    }

    if ( (*iter)->acSourceData_ != 0 )
    {
      (*iter)->acSourceData_->printOutParams ();
    }
  }

  os << std::endl;
  }
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
// ISRC Master functions:
//-----------------------------------------------------------------------------

Master::Master(
  const Configuration & configuration,
  const FactoryBlock &  factory_block,
  const SolverState &   solver_state,
  const DeviceOptions & device_options)
  : DeviceMaster<Traits>(configuration, factory_block, solver_state, device_options),
    HBSpecified_(factory_block.deviceManager_.getHBSpecified()),
    ACSpecified_(factory_block.deviceManager_.getACSpecified())
{}

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
    Instance & inst = *(*it);

    SourceData *dataPtr = inst.dcSourceData_; // by default assume the DC value.


    if ( (getSolverState().tranopFlag || (ACSpecified_ && !inst.DCSOURCETYPEgiven ) ) && getSolverState().locaEnabledFlag && inst.tranSourceData_ != 0)
    {

      if (inst.firstTimeload)
      {
        double  val = inst.tranSourceData_->returnSource();

        inst.setParam("DCV0", val,  true);

        inst.dcSourceData_->setParams (&inst.DCV0);

        inst.firstTimeload = false;
      }
    }

    if ((HBSpecified_ || (( getSolverState().tranopFlag || (ACSpecified_ && !inst.DCSOURCETYPEgiven ) )   && !getSolverState().locaEnabledFlag ) || getSolverState().transientFlag ) && inst.tranSourceData_ != 0 )


//    if ((HBSpecified_ || (getSolverState().tranopFlag && (!getSolverState().locaEnabledFlag ||  !inst.DCSOURCETYPEgiven ) ) || getSolverState().transientFlag || (ACSpecified_ && !inst.DCSOURCETYPEgiven ) ) && inst.tranSourceData_ != 0 )
    {
      dataPtr = inst.tranSourceData_;
    }

    double source = 0.0;
    if (dataPtr != 0)
    {
      source = dataPtr->returnSource();
    }

    bVec[inst.li_Pos] -= source;
    bVec[inst.li_Neg] += source;

    if( inst.loadLeadCurrent )
    {
      leadF[inst.li_branch_data] = source;
      junctionV[inst.li_branch_data] = solVec[inst.li_Pos] - solVec[inst.li_Neg];
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
  return true;
}


Device *Traits::factory(const Configuration &configuration, const FactoryBlock &factory_block)
{

  return new Master(configuration, factory_block, factory_block.solverState_, factory_block.deviceOptions_);
}

void
registerDevice(const DeviceCountMap& deviceMap, const std::set<int>& levelSet)
{
  if (deviceMap.empty() || (deviceMap.find("I")!=deviceMap.end()))
  {
    Config<Traits>::addConfiguration()
      .registerDevice("i", 1);
  }
}

} // namespace Resistor
} // namespace Device
} // namespace Xyce

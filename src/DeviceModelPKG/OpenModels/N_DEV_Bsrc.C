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
// Purpose        :
//
// Special Notes  :
//
// Creator        : Robert J Hoekstra, SNL, Parallel Computational Sciences
//
// Creation Date  : 06/05/01
//
//
//
//
//-------------------------------------------------------------------------
#include <Xyce_config.h>

// ---------- Standard Includes ----------
#include <N_UTL_Math.h>

// ----------   Xyce Includes   ----------
#include <N_DEV_Bsrc.h>
#include <N_DEV_DeviceOptions.h>
#include <N_DEV_ExternData.h>
#include <N_DEV_MatrixLoadData.h>
#include <N_DEV_SolverState.h>
#include <N_DEV_Message.h>
#include <N_ERH_ErrorMgr.h>

#include <N_LAS_Vector.h>
#include <N_LAS_Matrix.h>

#include <N_UTL_Expression.h>
#include <N_UTL_FeatureTest.h>

namespace Xyce {
namespace Device {


namespace Bsrc {


void Traits::loadInstanceParameters(ParametricData<Bsrc::Instance> &p)
{
  // Set up configuration constants:
  // Set up double precision variables:
  p.addPar("I", 0.0, &Bsrc::Instance::I)
    .setExpressionAccess(ParameterType::SOLN_DEP)
    .setUnit(U_AMP)
    .setDescription("Current for current source");

  p.addPar("V", 0.0, &Bsrc::Instance::V)
    .setExpressionAccess(ParameterType::SOLN_DEP)
    .setUnit(U_VOLT)
    .setDescription("Voltage for voltage source");

  p.addPar("TEMP", 0.0, &Bsrc::Instance::temp)
    .setUnit(U_DEGC)
    .setCategory(CAT_NONE)
    .setDescription("Device temperature");

  p.addPar ("SMOOTHBSRC",false, &Bsrc::Instance::newABM)
   .setUnit(U_LOGIC)
   .setCategory(CAT_NONE)
   .setDescription("smooth bsrc");

  p.addPar ("RCCONST", 1e-9, &Bsrc::Instance:: rc_const )
   .setUnit(U_SECOND)
   .setCategory(CAT_NONE)
   .setDescription("rc time constant");
}

void Traits::loadModelParameters(ParametricData<Bsrc::Model> &p)
{}


#define Xyce_NONPOINTER_MATRIX_LOAD 1

// Class Instance
//-----------------------------------------------------------------------------
// Function      : Instance::Instance
// Purpose       : "instance block" constructor
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 06/05/01
//-----------------------------------------------------------------------------
Instance::Instance(
  const Configuration & configuration,
  const InstanceBlock &         IB,
  Model &                       BMiter,
  const FactoryBlock &          factory_block)
  : DeviceInstance(IB, configuration.getInstanceParameters(), factory_block),
    model_(BMiter),
    Exp_ptr(0),
    expNumVars(0),
    expBaseVar(0),
    expNumDdt(0),
    expVal(0),
    IB(IB),
    isVSRC(false),
    scale(1.0),
    nlstep(-1),
    li_Pos(-1),
    li_Neg(-1),
    li_Bra(-1),
    li_branch_data(0),

    ABraEquPosNodeOffset(-1),
    ABraEquNegNodeOffset(-1),
    APosEquBraVarOffset(-1),
    ANegEquBraVarOffset(-1),
    rc_const (1e-9),
    APosEquPosNodeOffset(-1),
    ANegEquPosNodeOffset(-1),
    APosEquNegNodeOffset(-1),
    ANegEquNegNodeOffset(-1),
    newABM(0),
    fBraEquPosNodePtr(0),
    fBraEquNegNodePtr(0),
    fPosEquBraVarPtr(0),
    fNegEquBraVarPtr(0)
{
  numIntVars   = 1;
  numExtVars   = 2;
  numStateVars = 0;
  setNumBranchDataVars(0);             // by default don't allocate space in branch vectors
  numBranchDataVarsIfAllocated = 1;    // this is the space to allocate if lead current or power is needed.

  // Set params to constant default values:
  setDefaultParams ();

  // Set params according to instance line and constant defaults from metadata:
  setParams (IB.params);

  // Set any non-constant parameter defaults:
  if (given("I") && !given("V"))
  {
    isVSRC = false;
    // current source doesn't have current as part of the solution vector
    // so store it in the store vector
    setNumStoreVars(1);
  }
  else if (!given("I") && given("V"))
  {
    isVSRC = true;
  }
  else
  {
    UserError(*this) << "Must supply one of V= or I=";
  }



  if (!given("RCCONST"))
    rc_const = getDeviceOptions().rc_const;

  if (!given("SMOOTHBSRC"))
    newABM = getDeviceOptions().newABM;

  if (isVSRC)
  {

    if ( newABM )
      numIntVars = 0;
    else
      numIntVars = 1;

  }
  else
  {
    numIntVars = 0;
  }

  std::vector<Depend>::const_iterator d;
  std::vector<Depend>::const_iterator begin = getDependentParams().begin();
  std::vector<Depend>::const_iterator end = getDependentParams().end();

  for  (d = begin ; d != end ; ++d)
  {
    if (d->name == "I" || d->name == "V")
    {
      expNumVars = d->n_vars;
      expBaseVar = d->lo_var;
      Exp_ptr = d->expr;

      expNumDdt = Exp_ptr->getNumDdt();
      ddtVals.resize(expNumDdt);
      li_ddt.resize(expNumDdt);
      numStateVars += expNumDdt;

      if (expNumVars>0)
      {
        expVarDerivs.resize(expNumVars);
        myVarVals.resize(expNumVars);

        // this tells the device entity class NOT to call evaluateFunction on 
        // this expression, as that will be redundant with the evaluate call that 
        // is done from the Bsrc.  The exception is if the expression contains the 
        // ddt operator.  Then, it needs to call both  evaluateFunction and later 
        // evaluate for it to work properly.
        if (expNumDdt<=0)
        {
          dependentParamExcludeMap_[d->name] = 1;
        }
      }
      break;
    }
  }

  if( jacStamp.empty() )
  {
    if( isVSRC )
    {
      if ( newABM )
      {
        jacStamp.resize(2);
        jacStamp[0].resize(2 +expNumVars);
        jacStamp[1].resize(2 +expNumVars);
        jacStamp[0][0]=0;
        jacStamp[0][1] = 1;
        jacStamp[1][0] = 0;
        jacStamp[1][1] = 1;
        for( int i = 0; i < expNumVars; ++i )
        {
          jacStamp[0][i+2] = i+2;
          jacStamp[1][i+2] = i+2;
        }
      }
      else
      {
        jacStamp.resize(3);
        jacStamp[0].resize(1);
        jacStamp[0][0]=2;
        jacStamp[1].resize(1);
        jacStamp[1][0]=2;
        jacStamp[2].resize(2+expNumVars);
        jacStamp[2][0]=0;
        jacStamp[2][1]=1;
        for( int i = 0; i < expNumVars; ++i )
        {
          jacStamp[2][i+2] = i+3;         
        }
      }
    }
    else
    {
      jacStamp.resize( 2 );
      jacStamp[0].resize(expNumVars);
      jacStamp[1].resize(expNumVars);
      for( int i = 0; i < expNumVars; ++i )
      {
        jacStamp[0][i] = i+2;
        jacStamp[1][i] = i+2;
      }
    }
  }

  // Calculate any parameters specified as expressions:

  updateDependentParameters();

  // calculate dependent (ie computed) params and check for errors:

  processParams();

}

//-----------------------------------------------------------------------------
// Function      : Instance::processParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 07/30/03
//-----------------------------------------------------------------------------
bool Instance::processParams()
{
  updateTemperature(temp);

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::updateTemperature
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
bool Instance::updateTemperature(const double & temp_tmp)
{
  if (expNumVars == 0)
  {
    if (isVSRC)
    {
      expVal = V;
    }
    else
    {
      expVal = I;
    }
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::~Instance
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 06/05/01
//-----------------------------------------------------------------------------
Instance::~Instance ()
{
}

// Additional Declarations
//-----------------------------------------------------------------------------
// Function      : Instance::registerLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 6/21/02
//-----------------------------------------------------------------------------
void Instance::registerLIDs( const std::vector<int> & intLIDVecRef,
                             const std::vector<int> & extLIDVecRef )
{
  AssertLIDs(intLIDVecRef.size() == numIntVars);
  AssertLIDs(extLIDVecRef.size() == numExtVars);

  if (DEBUG_DEVICE)
  {
    Xyce::dout() << std::endl << section_divider << std::endl;
    Xyce::dout() << "  BsrcInstance::registerLIDs" << std::endl;
    Xyce::dout() << "  name = " << getName() << std::endl;
  }

  // copy over the global ID lists.
  intLIDVec = intLIDVecRef;
  extLIDVec = extLIDVecRef;

  // Now use these lists to obtain the indices into the
  // linear algebra entities.  This assumes an order.
  // For the matrix  indices, first do the rows.

  li_Pos = extLIDVec[0];
  li_Neg = extLIDVec[1];

  if (DEBUG_DEVICE)
  {
    Xyce::dout() << "  li_Pos = " << li_Pos << std::endl;
    Xyce::dout() << "  li_Neg = " << li_Neg << std::endl;
  }

  if( isVSRC && ! newABM )
  {
    li_Bra = intLIDVec[0];

    if (DEBUG_DEVICE)
    {
      Xyce::dout() << "  li_Bra = " << li_Bra << std::endl;
    }
  }     

  if (DEBUG_DEVICE)
  {
    Xyce::dout() << section_divider << std::endl;
  }
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
/// The B source device uses exactly one "branch data vector" element, where
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
  if (isVSRC)
    addInternalNode(symbol_table, li_Bra, getName(), "branch");

  if (loadLeadCurrent)
  {
    addBranchDataNode( symbol_table, li_branch_data, getName(), "BRANCH_D");
  }

}

//-----------------------------------------------------------------------------
// Function      : Instance::registerStoreLIDs
// Purpose       : One store var for device current if this is a current source
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical Systems Modeling
// Creation Date : 01/17/2013
//-----------------------------------------------------------------------------
void Instance::registerStoreLIDs(const std::vector<int> & stoLIDVecRef )
{
}

//-----------------------------------------------------------------------------
// Function      : Instance::registerStateLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 06/21/02
//-----------------------------------------------------------------------------
void Instance::registerStateLIDs( const std::vector<int> & staLIDVecRef )
{
  AssertLIDs(staLIDVecRef.size() == numStateVars);
  AssertLIDs(li_ddt.size() == expNumDdt);
  AssertLIDs(numStateVars == expNumDdt);

  for (int i=0 ; i<expNumDdt ; ++i)
  {
    li_ddt[i] = staLIDVecRef[i];
  }
}


//-----------------------------------------------------------------------------
// Function      : Instance::getDepSolnVars
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 06/06/01
//-----------------------------------------------------------------------------
const std::vector<std::string> & Instance::getDepSolnVars()
{
  return DeviceInstance::getDepSolnVars();
}


//-----------------------------------------------------------------------------
// Function      : Instance::jacobianStamp
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 9/2/02
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
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 9/2/02
//-----------------------------------------------------------------------------
void Instance::registerJacLIDs(
  const std::vector< std::vector<int> > & jacLIDVec)
{
  DeviceInstance::registerJacLIDs( jacLIDVec );
  if( isVSRC )
  {

    if ( newABM )
    {
      APosEquPosNodeOffset = jacLIDVec[0][0];
      APosEquNegNodeOffset = jacLIDVec[0][1];
      ANegEquPosNodeOffset = jacLIDVec[1][0];
      ANegEquNegNodeOffset = jacLIDVec[1][1];

      APosEquExpVarOffsets.resize( expNumVars );
      ANegEquExpVarOffsets.resize( expNumVars );

      for( int i = 0; i < expNumVars; ++i )
      {
        APosEquExpVarOffsets[i] = jacLIDVec[0][i+2];
        ANegEquExpVarOffsets[i] = jacLIDVec[1][i+2];
      }

    }
    else
    {
      APosEquBraVarOffset = jacLIDVec[0][0];
      ANegEquBraVarOffset = jacLIDVec[1][0];
      ABraEquPosNodeOffset = jacLIDVec[2][0];
      ABraEquNegNodeOffset = jacLIDVec[2][1];
      ABraEquExpVarOffsets.resize( expNumVars );
      for( int i = 0; i < expNumVars; ++i )
      {
        ABraEquExpVarOffsets[i] = jacLIDVec[2][i+2]; 
      }
    }
  }
  else
  {
    APosEquExpVarOffsets.resize( expNumVars );
    ANegEquExpVarOffsets.resize( expNumVars );
    for( int i = 0; i < expNumVars; ++i )
    {
      APosEquExpVarOffsets[i] = jacLIDVec[0][i];
      ANegEquExpVarOffsets[i] = jacLIDVec[1][i];
    }

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

  if( isVSRC )
  {
    fPosEquBraVarPtr = &(dFdx[li_Pos][APosEquBraVarOffset]);
    fNegEquBraVarPtr = &(dFdx[li_Neg][ANegEquBraVarOffset]);
    fBraEquPosNodePtr = &(dFdx[li_Bra][ABraEquPosNodeOffset]);
    fBraEquNegNodePtr = &(dFdx[li_Bra][ABraEquNegNodeOffset]);

    fBraEquExpVarPtrs.resize( expNumVars );
    for( int i = 0; i < expNumVars; ++i )
    {
      fBraEquExpVarPtrs[i] = &(dFdx[li_Bra][ ABraEquExpVarOffsets[i] ]);
    }
  }
  else
  {
    fPosEquExpVarPtrs.resize( expNumVars );
    fNegEquExpVarPtrs.resize( expNumVars );
    for( int i = 0; i < expNumVars; ++i )
    {
      fPosEquExpVarPtrs[i]  = &(dFdx[li_Pos][ APosEquExpVarOffsets[i] ]);
      fNegEquExpVarPtrs[i]  = &(dFdx[li_Neg][ ANegEquExpVarOffsets[i] ]);
    }
  }

#endif
}

//-----------------------------------------------------------------------------
// Function      : Instance::updateIntermediateVars
//
// Purpose       : Calls the expression handler to evaluate the expression
//                 and various derivatives.  These quantities are needed
//                 for the vector and matrix loads.
//
// Special Notes :
//
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 06/05/01
//-----------------------------------------------------------------------------
bool Instance::updateIntermediateVars ()
{
  if (expNumVars == 0)
  {
    if (isVSRC)
    {
      expVal = V;
    }
    else
    {
      expVal = I;
    }
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::updatePrimaryState
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 06/05/01
//-----------------------------------------------------------------------------
bool Instance::updatePrimaryState ()
{
  bool bsuccess=updateIntermediateVars ();

  // Get values of the arguments for ddt() calls in expression so that the derivatives
  // can be determined by the time integration class
  if (expNumDdt > 0)
  {
    double * staVec = extData.nextStaVectorRawPtr;

    Exp_ptr->getDdtVals (ddtVals);
    for (int i=0 ; i<expNumDdt ; ++i)
    {
      staVec[li_ddt[i]] = ddtVals[i];
    }
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::updateSecondaryState
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 06/05/01
//-----------------------------------------------------------------------------
bool Instance::updateSecondaryState ()
{
  // Get time derivatives from time integrator, and evaluate expression to get
  // derivatives with respect to independent quantities

  if (expNumDdt > 0)
  {
    double * staDerivVec = extData.nextStaDerivVectorRawPtr;

    for (int i=0 ; i<expNumDdt ; ++i)
    {
      ddtVals[i] = staDerivVec[li_ddt[i]];
    }
    Exp_ptr->setDdtDerivs(ddtVals);
  }
  // Evaluate Expression with corrected time derivative values
  if (expNumVars != 0)
  {
    Exp_ptr->evaluate( expVal, expVarDerivs);
  }

  // Test derivatives, if too big, zero out
  for (int i = 0; i < expNumVars; ++i)
  {
    double maxMag = 1.0e+10;
    if (expVarDerivs[i] > maxMag || expVarDerivs[i] < -maxMag)
    {
      static Report::MessageCode id;

      Report::UserWarning(id) << "In device " << getName() << ": Expression derivative for variable number " << i << " |" << expVarDerivs[i] << "| exceeds " << maxMag << ", value reduced";

      expVarDerivs[i] = (expVarDerivs[i] > 0) ? maxMag : -maxMag;
    }
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEFVector
//
// Purpose       : Loads the F-vector contributions for a single
//                 bsrc instance.
//
// Special Notes : See the special notes for loadDAEFVector.
//
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 04/27/04
//-----------------------------------------------------------------------------
bool Instance::loadDAEFVector ()
{
  double source(0.0), v_pos(0.0), v_neg(0.0), i_bra(0.0);
  double * solVec = extData.nextSolVectorRawPtr;
  double * fVec = extData.daeFVectorRawPtr;
  double * stoVec = extData.nextStoVectorRawPtr;

  source = expVal;

  //VSRC or ISRC
  if (isVSRC)
  {
    // get the value for v_pos, v_neg, i_bra
    v_pos = solVec[li_Pos];
    v_neg = solVec[li_Neg];
    i_bra = solVec[li_Bra];

    double c_tmp = i_bra;
    double v_tmp = (v_pos-v_neg-source);

    fVec[li_Pos] +=  c_tmp;
    fVec[li_Neg] += -c_tmp;
    fVec[li_Bra] +=  v_tmp;

    if( loadLeadCurrent )
    {
      double * leadF = extData.nextLeadCurrFCompRawPtr;
      leadF[li_branch_data] = c_tmp;
      double * junctionV = extData.nextJunctionVCompRawPtr;
      junctionV[li_branch_data] = v_tmp; 
    }

  }
  else
  {
    fVec[li_Pos] += source;
    fVec[li_Neg] += -source;

    if( loadLeadCurrent )
    {
      double * leadF = extData.nextLeadCurrFCompRawPtr;
      leadF[li_branch_data] = source;
      double * junctionV = extData.nextJunctionVCompRawPtr;
      junctionV[li_branch_data] = solVec[li_Pos] - solVec[li_Neg];
    }
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEdFdx ()
//
// Purpose       : Loads the F-vector contributions for a single
//                 bsrc instance.
//
// Special Notes :
//
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 04/27/04
//-----------------------------------------------------------------------------
bool Instance::loadDAEdFdx ()
{
  Linear::Matrix & dFdx = *(extData.dFdxMatrixPtr);

  double coef = 1.0;

  if( isVSRC )
  {
    dFdx[li_Pos][APosEquBraVarOffset]  += coef;
    dFdx[li_Neg][ANegEquBraVarOffset]  -= coef;
    dFdx[li_Bra][ABraEquPosNodeOffset] += coef;
    dFdx[li_Bra][ABraEquNegNodeOffset] -= coef;

    for( int i = 0; i < expNumVars; ++i )
    {
      dFdx[li_Bra][ABraEquExpVarOffsets[i]] -= expVarDerivs[i];
    }
  }
  else
  {
    if( expNumVars )
    {
      for( int i = 0; i < expNumVars; ++i )
      {
        dFdx[li_Pos][APosEquExpVarOffsets[i]] += expVarDerivs[i];
        dFdx[li_Neg][ANegEquExpVarOffsets[i]] -= expVarDerivs[i];
      }
    }
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::varTypes
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 02/17/02
//-----------------------------------------------------------------------------
void Instance::varTypes( std::vector<char> & varTypeVec )
{
  // if the B-Source has a Voltage-source formulation then it will contribute 
  // a branch current to the solution vector.
  if( isVSRC )
  {
    varTypeVec.resize(1);
    varTypeVec[0] = 'I';
  }
}

// Class Model

//-----------------------------------------------------------------------------
// Function      : Model::Model
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 06/05/01
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
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 06/05/01
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
// Function      : DeviceModel::processParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 06/03/02
//-----------------------------------------------------------------------------
bool Model::processParams()
{
  return true;
}

//-----------------------------------------------------------------------------
// Function      : DeviceModel::processInstanceParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 03/23/06
//-----------------------------------------------------------------------------
bool Model::processInstanceParams()
{
  return true;
}

// Additional Declarations

//-----------------------------------------------------------------------------
// Function      : Model::printOutInstances
// Purpose       : debugging tool.
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 6/05/01
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

// Bsrc Master functions:

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
  for (InstanceVector::const_iterator it = getInstanceBegin(); it != getInstanceEnd(); ++it)
  {
    Instance & bi =  *(*it);
    if (bi.expNumVars == 0)
    {
      if (bi.isVSRC)
      {
        bi.expVal = bi.V;
      }
      else
      {
        bi.expVal = bi.I;
      }
    }
    // Get values of the arguments for ddt() calls in expression so that the derivatives
    // can be determined by the time integration class
    if (bi.expNumDdt > 0)
    {
      bi.Exp_ptr->getDdtVals (bi.ddtVals);
      for (int j=0 ; j<bi.expNumDdt ; ++j)
      {
        staVec[bi.li_ddt[j]] = bi.ddtVals[j];
      }
    }
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Master::updateSecondaryState
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 11/26/08
//-----------------------------------------------------------------------------
bool Master::updateSecondaryState ( double * staDerivVec, double * stoVec )
{
  for (InstanceVector::const_iterator it = getInstanceBegin(); it != getInstanceEnd(); ++it)
  {
    Instance & bi =  *(*it);

    // Get time derivatives from time integrator, and evaluate expression to get
    // derivatives with respect to independent quantities

    if (bi.expNumDdt > 0)
    {
      for (int j=0 ; j<bi.expNumDdt ; ++j)
      {
        bi.ddtVals[j] = staDerivVec[bi.li_ddt[j]];
      }
      bi.Exp_ptr->setDdtDerivs(bi.ddtVals);
    }
    // Evaluate Expression with corrected time derivative values
    if (bi.expNumVars != 0)
    {
      bi.Exp_ptr->evaluate( bi.expVal, bi.expVarDerivs);
    }

    // Test derivatives, if too big, zero out
    for (int k = 0; k < bi.expNumVars; ++k)
    {
      double maxMag = 1.0e+10;
      if (bi.expVarDerivs[k] > maxMag || bi.expVarDerivs[k] < -maxMag)
      {
        static Report::MessageCode id;

        Report::UserWarning(id) << "In device " << bi.getName() << ": Expression derivative for variable number " << k << " |" << bi.expVarDerivs[k] << "| exceeds " << maxMag << ", value reduced";
        bi.expVarDerivs[k] = (bi.expVarDerivs[k] > 0) ? maxMag : -maxMag;
      }
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
    Instance & bi = *(*it);
    double v_pos(0.0), v_neg(0.0), i_bra(0.0);
    double source = bi.expVal;

    //VSRC or ISRC
    if (bi.isVSRC)
    {
      v_pos = solVec[bi.li_Pos];
      v_neg = solVec[bi.li_Neg];

      if ( bi.newABM)
      {
        double v_tmp1 = (v_pos-v_neg);
        double c_tmp1 = v_tmp1/1e6;

        double v_tmp = (v_pos-v_neg-source);

        double q_tmp = (v_pos - v_neg)*( bi.rc_const/1e-3);

        double c_tmp = v_tmp/1e-3;

        fVec[bi.li_Pos] +=  c_tmp +  c_tmp1;
        fVec[bi.li_Neg] += -c_tmp - c_tmp1;

        qVec[bi.li_Pos] +=  q_tmp;
        qVec[bi.li_Neg] += -q_tmp;

        // get the value for v_pos, v_neg, i_bra
      }
      else
      {
        i_bra = solVec[bi.li_Bra];

        double c_tmp = i_bra;
        double v_tmp = (v_pos-v_neg-source);

        fVec[bi.li_Pos] +=  c_tmp;
        fVec[bi.li_Neg] += -c_tmp;
        fVec[bi.li_Bra] +=  v_tmp;

        if( bi.loadLeadCurrent )
        {
          leadF[bi.li_branch_data] = c_tmp;
          junctionV[bi.li_branch_data] = v_pos-v_neg; 
        }
      }
    }
    else
    {
      fVec[bi.li_Pos] += source;
      fVec[bi.li_Neg] += -source;
      if( bi.loadLeadCurrent )
      {
        leadF[bi.li_branch_data] = source;
        junctionV[bi.li_branch_data] = solVec[bi.li_Pos] - solVec[bi.li_Neg];
      }
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
    Instance & bi = *(*it);
    double coef = 1.0;

    if( bi.isVSRC )
    {
#ifndef Xyce_NONPOINTER_MATRIX_LOAD
      *bi.fPosEquBraVarPtr  += coef;
      *bi.fNegEquBraVarPtr  -= coef;
      *bi.fBraEquPosNodePtr += coef;
      *bi.fBraEquNegNodePtr -= coef;

      for( int j = 0; j < bi.expNumVars; ++j )
      {
        *bi.fBraEquExpVarPtrs[j] -= bi.expVarDerivs[j];
      }
#else

//      coef = 1/1e-2 + 1/1e6;

      if ( bi.newABM == 1)
      {
        coef = 1/1e-3 + 1/1e6;

        dFdx[bi.li_Pos][bi.APosEquPosNodeOffset]  += coef;
        dFdx[bi.li_Pos][bi.APosEquNegNodeOffset]  -= coef;
        dFdx[bi.li_Neg][bi.ANegEquPosNodeOffset]  -= coef;
        dFdx[bi.li_Neg][bi.ANegEquNegNodeOffset]  += coef;

        for( int j = 0; j < bi.expNumVars; ++j )
        {
          dFdx[bi.li_Pos][bi.APosEquExpVarOffsets[j]]  +=  - bi.expVarDerivs[j]/1e-3;
          dFdx[bi.li_Neg][bi.ANegEquExpVarOffsets[j]]  += bi.expVarDerivs[j]/1e-3;
        }

        double c = bi.rc_const/ 1e-3;
        dQdx[bi.li_Pos][bi.APosEquPosNodeOffset]  += c ;
        dQdx[bi.li_Pos][bi.APosEquNegNodeOffset]  -= c;
        dQdx[bi.li_Neg][bi.ANegEquPosNodeOffset]  -= c;
        dQdx[bi.li_Neg][bi.ANegEquNegNodeOffset]  += c;
      }
      else
      {
        dFdx[bi.li_Pos][bi.APosEquBraVarOffset]  += coef;
        dFdx[bi.li_Neg][bi.ANegEquBraVarOffset]  -= coef;
        dFdx[bi.li_Bra][bi.ABraEquPosNodeOffset] += coef;
        dFdx[bi.li_Bra][bi.ABraEquNegNodeOffset] -= coef;

        for( int j = 0; j < bi.expNumVars; ++j )
      {
        dFdx[bi.li_Bra][bi.ABraEquExpVarOffsets[j]] -= bi.expVarDerivs[j];
      }


      }

#endif
    }
    else
    {
      if( bi.expNumVars )
      {
#ifndef Xyce_NONPOINTER_MATRIX_LOAD
        for( int k = 0; k < bi.expNumVars; ++k )
        {
          *bi.fPosEquExpVarPtrs[k] += bi.expVarDerivs[k];
          *bi.fNegEquExpVarPtrs[k] -= bi.expVarDerivs[k];
        }
#else
        for( int j = 0; j < bi.expNumVars; ++j )
        {
          dFdx[bi.li_Pos][bi.APosEquExpVarOffsets[j]] += bi.expVarDerivs[j];
          dFdx[bi.li_Neg][bi.ANegEquExpVarOffsets[j]] -= bi.expVarDerivs[j];
        }
#endif
      }
    }
  }
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadErrorWeightMask
//
// Purpose       : Loads the zero elements of the device mask
//
// Special Notes : elements of the error vector associated with zero
//                 elements of the mask will not be included in weighted
//                 norms by the time integrator.
//
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 01/19/07
//-----------------------------------------------------------------------------
/*void Instance::loadErrorWeightMask ()
{
  Linear::Vector * maskVectorPtr = extData.deviceErrorWeightMask_;

//  (*maskVectorPtr)[li_Pos] = 0.0;

//  (*maskVectorPtr)[li_Neg] = 0.0;

}   */

Device *Traits::factory(const Configuration &configuration, const FactoryBlock &factory_block)
{
  return new Master(configuration, factory_block, factory_block.solverState_, factory_block.deviceOptions_);
}

void
registerDevice(const DeviceCountMap& deviceMap, const std::set<int>& levelSet)
{
  static bool initialized = false;

  if (!initialized && (deviceMap.empty() 
       || (deviceMap.find("B")!=deviceMap.end())
       || (deviceMap.find("F")!=deviceMap.end())
       || (deviceMap.find("H")!=deviceMap.end())))
  {
    initialized = true;

    Config<Traits>::addConfiguration()
      .registerDevice("b", 1)
      .registerDevice("f", 1)
      .registerDevice("h", 1);
  }
}

} // namespace Bsrc
} // namespace Device
} // namespace Xyce

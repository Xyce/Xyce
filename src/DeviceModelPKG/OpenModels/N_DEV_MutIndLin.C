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
#include <N_UTL_Math.h>
#include <algorithm>
#include <set>

// ----------   Xyce Includes   ----------
#include <N_DEV_DeviceOptions.h>
#include <N_DEV_DeviceMaster.h>
#include <N_DEV_ExternData.h>
#include <N_DEV_MatrixLoadData.h>
#include <N_DEV_MutIndLin.h>
#include <N_DEV_SolverState.h>
#include <N_DEV_Message.h>
#include <N_ERH_ErrorMgr.h>

#include <N_LAS_Vector.h>
#include <N_LAS_Matrix.h>

#include <N_UTL_Expression.h>
#include <N_UTL_FeatureTest.h>
#include <N_UTL_HspiceBools.h>

namespace Xyce {
namespace Device {


//-----------------------------------------------------------------------------
// Function      : InductorInstanceData::InductorInstanceData
// Purpose       :
// Special Notes :
//
// Need to have a constructor for this data-only class
// or it is treated as in-line by some compilers (gcc 3.3 on the mac)
// Trying to inline this could cause problems
//
// Scope         : public
// Creator       : Rich Schiek
// Creation Date :
//-----------------------------------------------------------------------------
InductorInstanceData::InductorInstanceData():
    name(""),    // name of inductor
    L(0.0),
    IC(0.0),
    ICGiven(false),
    baseL(0.0),
    li_Pos(-1),
    li_Neg(-1),
    li_Branch(-1),
    APosEquBraVarOffset(-1),
    ANegEquBraVarOffset(-1),
    ABraEquPosNodeOffset(-1),
    ABraEquNegNodeOffset(-1),
    ABraEquBraVarOffset(-1),
    magOffset(-1),
    vPosOffset(-1),
    vNegOffset(-1)
#ifndef Xyce_NONPOINTER_MATRIX_LOAD
    ,
    // Pointers for dFdx entries:
    f_PosEquBraVarPtr(0),
    f_NegEquBraVarPtr(0),
    f_BraEquPosNodePtr(0),
    f_BraEquNegNodePtr(0),
    f_BraEquBraVarPtr(0),
    // offsets only needed in nonlinear application
    f_magPtr(0),
    f_vPosPtr(0),
    f_vNegPtr(0),
    q_PosEquBraVarPtr(0),
    q_NegEquBraVarPtr(0),
    q_BraEquPosNodePtr(0),
    q_BraEquNegNodePtr(0),
    q_BraEquBraVarPtr(0),
    // offsets only needed in nonlinear application
    q_magPtr(0),
    q_vPosPtr(0),
    q_vNegPtr(0)
#endif
{}

namespace MutIndLin {

void Traits::loadInstanceParameters(ParametricData<MutIndLin::Instance> &p)
{
  p.addPar ("COUP_VAL",1.0,&MutIndLin::Instance::mutualCup)
   .setGivenMember(&MutIndLin::Instance::mutualCupGiven)
   .setUnit(U_NONE)
   .setCategory(CAT_NONE)
   .setDescription("Coupling value");

  p.addPar ("COUPLEDMutIndLin",std::vector<std::string>(),&MutIndLin::Instance::inductorNames)
   .setUnit(U_NONE)
   .setCategory(CAT_NONE)
   .setDescription("");

  p.addPar ("COUPLEDINDUCTANCE",std::vector<double>(),&MutIndLin::Instance::inductorInductances)
   .setUnit(U_NONE)
   .setCategory(CAT_NONE)
   .setDescription("");

  p.addPar ("NODE1",std::vector<std::string>(),&MutIndLin::Instance::inductorsNode1)
   .setUnit(U_NONE)
   .setCategory(CAT_NONE)
   .setDescription("");

  p.addPar ("NODE2",std::vector<std::string>(),&MutIndLin::Instance::inductorsNode2)
   .setUnit(U_NONE)
   .setCategory(CAT_NONE)
   .setDescription("");

  p.addPar ("COUPLING",std::vector<double>(),&MutIndLin::Instance::couplingCoefficient)
   .setExpressionAccess(ParameterType::SOLN_DEP)
   .setUnit(U_NONE)
   .setCategory(CAT_NONE)
   .setDescription("Coupling coefficient");

  p.addPar ("COUPLEDINDUCTOR",std::vector<std::string>(),&MutIndLin::Instance::couplingInductor)
   .setUnit(U_NONE)
   .setCategory(CAT_NONE)
   .setDescription("");
   
  p.addPar ("IC",std::vector<double>(),&MutIndLin::Instance::initialCondition)
   .setUnit(U_AMP)
   .setCategory(CAT_NONE)
   .setDescription("Initial current through the inductor.");
}

void Traits::loadModelParameters(ParametricData<MutIndLin::Model> &p)
{
  p.addPar ("TNOM",27.0,&MutIndLin::Model::tnom)
      .setUnit(U_DEGC)
      .setCategory(CAT_MATERIAL)
      .setDescription("Reference temperature");

  p.addPar ("TC1",0.0,&MutIndLin::Model::tempCoeff1)
      .setUnit(U_NONE)
      .setCategory(CAT_MATERIAL)
      .setDescription("First order temperature coeff.");

  p.addPar ("TC2",0.0,&MutIndLin::Model::tempCoeff2)
      .setUnit(U_NONE)
      .setCategory(CAT_MATERIAL)
      .setDescription("Second order temperature coeff.");
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
    numInductors(0),
    mutualCup(0.0),
    mutualCupGiven(false),
    temp(getDeviceOptions().temp.getImmutableValue<double>()),
    tempGiven(false),
    scalingRHS(1.0)
{
  scalingRHS = 1.0;

  // set some default values.  May be changed by processParams
  numExtVars   = 2;
  numIntVars   = 1;
  numStateVars = 0;

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
#ifndef Xyce_NONPOINTER_MATRIX_LOAD
    inductorData->f_inductorCurrentPtrs.resize( inductorNames.size() );
    inductorData->q_inductorCurrentPtrs.resize( inductorNames.size() );
#endif
    inductorData->depVarPairs.clear();
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
  for(int i=0; i<numInductors; i++)
  {
    devConMap[i] = devConMap[i+1] = (i+1);
  }

  // now assemble the mutual coupling coefficient matrix
  // assume no coupling except to self, so just ones on dialgonal.
  mutualCouplingCoef.resize( numInductors );
  mutualCouplingCoefDerivs.resize( numInductors );
  for( int i=0; i<numInductors; ++i )
  {
    mutualCouplingCoef[i].resize( numInductors );
    fill( mutualCouplingCoef[i].begin(), mutualCouplingCoef[i].end(), 0.0);
    mutualCouplingCoef[i][i] = 1.0;

    mutualCouplingCoefDerivs[i].resize( numInductors );
    fill( mutualCouplingCoefDerivs[i].begin(),
          mutualCouplingCoefDerivs[i].end(), 0.0 );
  }

  // take the vector couplingInductor in pairs to figure out
  // the coupling value between inductors then insert it into
  // the couplingCoefficient matrix.
  // must take sequential pairs to avoid missing some when there are an odd number of inductors.
  // for example L1 L2 L3 has pairs L1-L2, L1-L3, L2-L3 
  // however, depending on how the netlist was created, the couplingInductor list could be
  // "L1" "L2" "L1" "L3" ... even pairs of numbers
  // or
  // "L1" "L2" "L3" "L4" ... potentially even or odd but inductor names do not repeat.
  //
  bool foundDuplicateInductorName=false;
  if (couplingInductor.size() % 2 == 0)
  {
    for( int i=0; i<couplingInductor.size(); ++i)
    {
      for( int j=i+1; j<couplingInductor.size(); ++j)
      {
        if( (couplingInductor[i] == couplingInductor[j] ) )
        {
          foundDuplicateInductorName=true;
          break; 
        }
      }
    }
  }
  
  if( foundDuplicateInductorName )
  {
    // in this case the list of coupledInductors is a list of pairs as in
    // "L1" "L2" "L1" "L3" etc.  treat is as such.

    // take the vector couplingInductor in pairs to figure out
    // the coupling value between inductors then insert it into
    // the couplingCoefficient matrix.
    indexPairs.resize(couplingInductor.size()/2);
    for( int i=0; i<(couplingInductor.size()-1) ; i+=2 )
    {
      // from inductor names, derive indices.
      int indexInductorOne = -1;
      int indexInductorTwo = -1;
      for( int j=0; j<inductorNames.size(); ++j)
      {
        if( couplingInductor[i] == inductorNames[j] )
        {
          indexInductorOne = j;
        }
        if( couplingInductor[i+1] == inductorNames[j] )
        {
          indexInductorTwo = j;
        }
        if( indexInductorOne != -1 && indexInductorTwo != -1 )
        {
          // stop searching if we found the answers.
          break;
        }
      }
      // Save the indices for later
      indexPairs[i/2].first = indexInductorOne;
      indexPairs[i/2].second = indexInductorTwo;

      // with the indices at hand, assign the value.
      // note that couplingCoefficient can be of length 1 if all coefficients were the
      // same.  Catch that case here.
      if( (i/2) < couplingCoefficient.size() )
      {
        mutualCouplingCoef[ indexInductorOne ][ indexInductorTwo ] =
          mutualCouplingCoef[ indexInductorTwo ][ indexInductorOne ] = couplingCoefficient[ (i / 2) ];
      }
      else
      {
        mutualCouplingCoef[ indexInductorOne ][ indexInductorTwo ] =
          mutualCouplingCoef[ indexInductorTwo ][ indexInductorOne ] = couplingCoefficient[ 0 ];
      }
    }

  }
  else
  {
    // in this case the list of coupledInductors is a unique list as in 
    // "L1" "L2" "L3" etc.  treat it as such
    int numUniquePairs = couplingInductor.size()*(couplingInductor.size()-1)/2;
    indexPairs.resize(numUniquePairs);
    int indexPairCounter=0;
    for( int i=0; i<couplingInductor.size() ; i++ )
    {
      for( int j=(i+1); j<couplingInductor.size() ; j++ )
      {
        // from inductor names, derive indices.
        int indexInductorOne = -1;
        int indexInductorTwo = -1;
        for( int k=0; k<inductorNames.size(); ++k)
        {
          if( couplingInductor[i] == inductorNames[k] )
          {
            indexInductorOne = k;
          }
          if( couplingInductor[j] == inductorNames[k] )
          {
            indexInductorTwo = k;
          }
          if( indexInductorOne != -1 && indexInductorTwo != -1 )
          {
            // stop searching if we found the answers.
            break;
          }
        }
        // Save the indices for later
        indexPairs[indexPairCounter].first = indexInductorOne;
        indexPairs[indexPairCounter].second = indexInductorTwo;
        if( indexPairCounter < (numUniquePairs-1) )
          ++indexPairCounter;

        // with the indices at hand, assign the value.
        // note that couplingCoefficient can be of length 1 if all coefficients were the
        // same.  Catch that case here.
        if( (i/2) < couplingCoefficient.size() )
        {
          mutualCouplingCoef[ indexInductorOne ][ indexInductorTwo ] =
            mutualCouplingCoef[ indexInductorTwo ][ indexInductorOne ] = couplingCoefficient[ (i / 2) ];
        }
        else
        {
          mutualCouplingCoef[ indexInductorOne ][ indexInductorTwo ] =
            mutualCouplingCoef[ indexInductorTwo ][ indexInductorOne ] = couplingCoefficient[ 0 ];
        }
      }
    }
  }
  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
  {
    Xyce::dout() << "mutualCouplingCoef:  " << std::endl;
    for( int i=0; i<numInductors; ++i )
    {
      for( int j=0; j<numInductors; ++j )
      {
        Xyce::dout() << mutualCouplingCoef[i][j] << "  ";
      }
      Xyce::dout() << std::endl;
    }
    Xyce::dout() << "index pairs = ";
    for( int i=0; i<indexPairs.size() ; i++ )
    {
      Xyce::dout() << "( " << indexPairs[i].first << ", " << indexPairs[i].second << ") ";
    }
    Xyce::dout() << std::endl;
  }

  inductorCurrents.resize( numInductors );
  dIdt.resize( numInductors );

  inductanceVals.resize( numInductors );
  LOI.resize( numInductors );
  LO.resize( numInductors );
  for( int i=0; i<numInductors; ++i)
  {
    LO[i].resize( numInductors );
  }

  updateInductanceMatrix();

  // Calculate any parameters specified as expressions:
  updateDependentParameters();

  // calculate dependent (ie computed) params and check for errors:
  processParams ();

  // update internal/external/state variable counts
  numExtVars = 2*numInductors;
  numIntVars = numInductors;
  numStateVars = 0;

  numCoupDepVars.resize(couplingCoefficient.size(),0);
  expPtrs.resize(couplingCoefficient.size(),0); // null pointers
  couplingCoefficientVarDerivs.resize(couplingCoefficient.size());

  // set up the jacobian stamp
  // for an individual inductor the stamp would be:
  //
  //          V1   V2   Ib
  //  kcl1               1
  //  kcl2              -1
  //  branch  1    -1   L/dt
  //
  //  for a collection of these, the internal variable, branch equations,
  //  must be at the end of a given stamp row.
  //
  //  So for N inductors the samp is:
  //
  //           V1  V2  V3  V4 ... V2N  I1  I2  ... IN
  //  kcl1                              1
  //  kcl2                             -1
  //  kcl3                                  1
  //  kcl4                                 -1
  //  branch1  1   -1                 L/dt  c  ... c
  //  branch2           1  -1          c  L/dt ... c
  //
  // where "c" is an induced current change.

  jacStamp.resize( 3 * numInductors );
  for( int i=0; i< numInductors; ++i )
  {
    jacStamp[2*i].resize(1);          // Vpos row
    jacStamp[2*i+1].resize(1);        // Vneg row
    jacStamp[2*numInductors + i].resize(numInductors + 2); // Ibranch row

    jacStamp[2*i  ][0] = 2*numInductors + i;    // vpos-ibranch
    jacStamp[2*i+1][0] = 2*numInductors + i;    // vneg-ibranch
    jacStamp[2*numInductors + i][0] = 2*i;      // ibranch-vpos
    jacStamp[2*numInductors + i][1] = 2*i + 1;  // ibranch-vneg
    for( int j=0; j<numInductors; ++j )
    {
      jacStamp[2*numInductors + i][j+2] = 2*numInductors + j;
                                               // ibranch-ibranch[j]
    }
  }

  // Now we process our dependent variables.
  // We need to bump up the jacstamp for all the dependencies on
  // the coupling vars that are expressions.  We also need to keep track of
  // how many variables each inductor depends on.
  std::vector<Depend>::const_iterator d;
  std::vector<Depend>::const_iterator begin=getDependentParams().begin();
  std::vector<Depend>::const_iterator end=getDependentParams().end();

  for (d=begin; d != end; ++d)
  {

    if (d->name == "COUPLING" && d->vectorIndex != -1)
    {
      // sanity check, they should all have the first two conditions,
      // and the next means we actually have work to do
      if (d->n_vars >0)
      {
        // Keep track of the number of variables this coefficient depends on:
        numCoupDepVars[d->vectorIndex] = d->n_vars;

        // Save the pointer to the expresison, we need to be able to evaluate it.
        expPtrs[d->vectorIndex] = d->expr;

        // These are the two inductors that this coefficient couples
        int indexInductorOne = indexPairs[d->vectorIndex].first;
        int indexInductorTwo = indexPairs[d->vectorIndex].second;

        // The coupling coefficient introduces terms to both equations
        int inductorOneBranchSize = jacStamp[2*numInductors+
                                             indexInductorOne].size();
        int inductorTwoBranchSize = jacStamp[2*numInductors+
                                             indexInductorTwo].size();
        jacStamp[2*numInductors+indexInductorOne].resize(inductorOneBranchSize+
                                                         d->n_vars);
        jacStamp[2*numInductors+indexInductorTwo].resize(inductorTwoBranchSize+
                                                         d->n_vars);
        for (int i=0; i<d->n_vars; i++)
        {
          jacStamp[2*numInductors+indexInductorOne][inductorOneBranchSize+i]
            = d->lo_var+i+3*numInductors;
          jacStamp[2*numInductors+indexInductorTwo][inductorTwoBranchSize+i]
            = d->lo_var+i+3*numInductors;

          // We need to be able to determine which of our dependent var
          // indices map onto which variable of what coupling coefficient
          // when we go to do the jacobian loads.  We do this so we don't
          // need to keep all kinds of two-dimensional arrays like
          // mutualCoupling that are just copies of the one-dimensonal
          // arrays like couplingCoefficient.
          instanceData[indexInductorOne]->depVarPairs.push_back(
             std::pair<int,int>(d->vectorIndex,i));
          instanceData[indexInductorTwo]->depVarPairs.push_back(
             std::pair<int,int>(d->vectorIndex,i));
        }
      }
    }
    else if (d->name == "COUPLEDINDUCTANCE")
    {
      // this is fine --- Inductors are allowed to be dependent params on
      // global params, and we might get here
    }
    else
    {
      Report::UserError() <<"Error in mutual inductor constructor for " << getName().getEncodedName()
                          << ": Parameter " << d->name << " value " << d->expr->get_expression() << " fails sanity check.";
    }
  }


  // We're now done processing all the dependent variables, so we know
  // everything we need to know about expressions and derivative numbers.
  for (int i=0; i<couplingCoefficient.size(); ++i)
  {
    couplingCoefficientVarDerivs[i].resize(numCoupDepVars[i]);
  }

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
  {
    Xyce::dout() << "Instance::Instance-----------------" << std::endl;
    Xyce::dout() << "numExtVars = " << numExtVars << std::endl
      << "numIntVars = " << numIntVars << std::endl
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
  std::vector< InductorInstanceData* >::iterator currentInductor = instanceData.begin();
  std::vector< InductorInstanceData* >::iterator endInductor = instanceData.end();
  for ( ; currentInductor != endInductor ; ++currentInductor)
  {
    delete *currentInductor;
  }
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
    ++currentInductor;
  }

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
  {
    Xyce::dout() << "Instance::registerLIDs----------------------------" << std::endl;
    currentInductor = instanceData.begin();
    i=0;
    while( currentInductor != endInductor )
    {
      Xyce::dout() << "Inductor [ " << i++ << " ] "
           << "   li_Pos = " << (*currentInductor)->li_Pos
           << "   li_Neg = " << (*currentInductor)->li_Neg
           << "   li_Branch = " << (*currentInductor)->li_Branch << std::endl;
      ++currentInductor;
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
  std::string baseName = getSubcircuitName(getName());

  for (int i = 0; i < instanceData.size(); ++i) {
    addInternalNode(symbol_table, intLIDVec[i], getName(), instanceData[i]->name + "_branch");
    std::string branchInductorName = baseName; 
    if( branchInductorName != "" )
      branchInductorName += Xyce::Util::separator;
    branchInductorName += instanceData[i]->name;
    InstanceName bInductorIName = InstanceName( branchInductorName );
    std::string encodedName = spiceInternalName( bInductorIName, "branch");
    addInternalNode(symbol_table, intLIDVec[i], encodedName);
    if (loadLeadCurrent)
    {
      addBranchDataNode(symbol_table,instanceData[i]->li_branch_data, bInductorIName, "BRANCH_D");
    }
  }

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
}

//----------------------------------------------------------------------------- 
// Function      : Instance::registerBranchDataLIDs 
// Purpose       : This allows P() and W() to work for the component inductors
//               : of a mutual inductor.
// Special Notes : 
// Scope         : public 
// Creator       : Pete Sholander, SNL 
// Creation Date : 3/21/17 
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
  }
  std::vector< InductorInstanceData* >::iterator currentInductor = instanceData.begin();
  std::vector< InductorInstanceData* >::iterator endInductor = instanceData.end();
  int i = 0;
  while( currentInductor != endInductor )
  {
    (*currentInductor)->APosEquBraVarOffset  = jacLIDVec[ 2*i     ][ 0 ];
    (*currentInductor)->ANegEquBraVarOffset  = jacLIDVec[ 2*i + 1 ][ 0 ];
    (*currentInductor)->ABraEquPosNodeOffset = jacLIDVec[ 2*numInductors + i ][ 0 ];
    (*currentInductor)->ABraEquNegNodeOffset = jacLIDVec[ 2*numInductors + i ][ 1 ];
    for( int j=0; j<numInductors; ++j )
    {
      if( i == j )
      {
        (*currentInductor)->ABraEquBraVarOffset  = jacLIDVec[ 2*numInductors + i ][ j + 2 ];
      }
      (*currentInductor)->inductorCurrentOffsets[ j ] = jacLIDVec[ 2*numInductors + i ][ j + 2 ];
    }
    // Now do the parts that are for the dependent variables
    int numdepvars=(*currentInductor)->depVarPairs.size();
    (*currentInductor)->ABraEquDepVarOffsets.resize(numdepvars);
    for (int j = 0; j < numdepvars; ++j)
    {
      (*currentInductor)->ABraEquDepVarOffsets[j] =
        jacLIDVec[2*numInductors+i][numInductors+2+j];
    }
    ++currentInductor;
    ++i;
  }

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
  {
    Xyce::dout() << "Instance::registerJacLIDs--------------------------" << std::endl;
    currentInductor = instanceData.begin();
    i=0;
    while( currentInductor != endInductor )
    {
      Xyce::dout() << "Inductor [ " << i << " ] " << (*currentInductor)->name
           << "   APosEquBraVarOffset = " << (*currentInductor)->APosEquBraVarOffset
           << "   ANegEquBraVarOffset = " << (*currentInductor)->ANegEquBraVarOffset
           << "   ABraEquPosNodeOffset = " << (*currentInductor)->ABraEquPosNodeOffset
           << "   ABraEquNegNodeOffset = " << (*currentInductor)->ABraEquNegNodeOffset
           << "   ABraEquBraVarOffset = " << (*currentInductor)->ABraEquBraVarOffset << std::endl;
      Xyce::dout() << "\tInductor branch offsets = { ";
      for( int j=0; j<numInductors ; ++j )
      {
        Xyce::dout() << (*currentInductor)->inductorCurrentOffsets[ j ] << ", ";
      }
      Xyce::dout() << "} " << std::endl;
      ++i;
      ++currentInductor;
    }
  }
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
  Linear::Matrix & dQdx = *(extData.dQdxMatrixPtr);

  std::vector< InductorInstanceData* >::iterator currentInductor = instanceData.begin();
  std::vector< InductorInstanceData* >::iterator endInductor = instanceData.end();
  int i = 0;
  while( currentInductor != endInductor )
  {
    // dFdx pointers:
    (*currentInductor)->f_PosEquBraVarPtr  =
        &(dFdx[((*currentInductor)->li_Pos)] [((*currentInductor)->APosEquBraVarOffset)] );

    (*currentInductor)->f_NegEquBraVarPtr =
        &(dFdx[((*currentInductor)->li_Neg)]   [((*currentInductor)->ANegEquBraVarOffset)] );

    (*currentInductor)->f_BraEquPosNodePtr =
        &(dFdx[((*currentInductor)->li_Branch)][((*currentInductor)->ABraEquPosNodeOffset)] );

    (*currentInductor)->f_BraEquNegNodePtr =
        &(dFdx[((*currentInductor)->li_Branch)][((*currentInductor)->ABraEquNegNodeOffset)] );
        
    (*currentInductor)->f_BraEquBraVarPtr =
        &(dFdx[((*currentInductor)->li_Branch)][((*currentInductor)->ABraEquBraVarOffset)] );

    for( int j=0; j<numInductors; ++j )
    {
      if( i == j )
      {
        //(*currentInductor)->ABraEquBraVarOffset  =
        // Is this only used for JMat?  (ie old DAE)?

      }
      (*currentInductor)->f_inductorCurrentPtrs[ j ] =
            &(dFdx[((*currentInductor)->li_Branch)][(*currentInductor)->inductorCurrentOffsets[j]] );
    }
    // Now do the parts that are for the dependent variables
    int numdepvars=(*currentInductor)->depVarPairs.size();
    (*currentInductor)->f_BraEquDepVarPtrs.resize(numdepvars);

    for (int j = 0; j < numdepvars; ++j)
    {
      (*currentInductor)->f_BraEquDepVarPtrs[j] =
            &(dFdx[((*currentInductor)->li_Branch)][(*currentInductor)->ABraEquDepVarOffsets[j]] );

    }

    // dQdx pointers:
    (*currentInductor)->q_PosEquBraVarPtr  =
        &(dQdx[((*currentInductor)->li_Pos)] [((*currentInductor)->APosEquBraVarOffset)] );

    (*currentInductor)->q_NegEquBraVarPtr =
        &(dQdx[((*currentInductor)->li_Neg)]   [((*currentInductor)->ANegEquBraVarOffset)] );

    (*currentInductor)->q_BraEquPosNodePtr =
        &(dQdx[((*currentInductor)->li_Branch)][((*currentInductor)->ABraEquPosNodeOffset)] );

    (*currentInductor)->q_BraEquNegNodePtr =
        &(dQdx[((*currentInductor)->li_Branch)][((*currentInductor)->ABraEquNegNodeOffset)] );
        
    (*currentInductor)->q_BraEquBraVarPtr =
        &(dQdx[((*currentInductor)->li_Branch)][((*currentInductor)->ABraEquBraVarOffset)] );

    for( int j=0; j<numInductors; ++j )
    {
      if( i == j )
      {
        //(*currentInductor)->ABraEquBraVarOffset  =
        // Is this only used for JMat?  (ie old DAE)?

      }
      (*currentInductor)->q_inductorCurrentPtrs[ j ] =
            &(dQdx[((*currentInductor)->li_Branch)][(*currentInductor)->inductorCurrentOffsets[j]] );
    }
    // Now do the parts that are for the dependent variables
    numdepvars=(*currentInductor)->depVarPairs.size();
    (*currentInductor)->q_BraEquDepVarPtrs.resize(numdepvars);

    for (int j = 0; j < numdepvars; ++j)
    {
      (*currentInductor)->q_BraEquDepVarPtrs[j] =
            &(dQdx[((*currentInductor)->li_Branch)][(*currentInductor)->ABraEquDepVarOffsets[j]] );

    }

    ++currentInductor;
    ++i;
  }



#endif
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
  int i, j;

  // take the vector couplingInductor in pairs to figure out
  // the coupling value between inductors then insert it into
  // the couplingCoefficient matrix.
  j = indexPairs.size();
  for( i=0; i<j ; ++i )
  {
    // note that couplingCoefficient can be of length 1 if all coefficients were the
    // same.  Catch that case here.
    if( i < couplingCoefficient.size() )
    {
      mutualCouplingCoef[ indexPairs[i].first ][ indexPairs[i].second ] =
        mutualCouplingCoef[ indexPairs[i].second ][ indexPairs[i].first ] = couplingCoefficient[i];
    }
    else
    {
      mutualCouplingCoef[ indexPairs[i].first ][ indexPairs[i].second ] =
        mutualCouplingCoef[ indexPairs[i].second ][ indexPairs[i].first ] = couplingCoefficient[0];
    }
  }
  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS)  && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << "processParams: mutualCouplingCoef:  " << std::endl;
    for( int i=0; i<numInductors; ++i )
    {
      for( int j=0; j<numInductors; ++j )
      {
        Xyce::dout() << mutualCouplingCoef[i][j] << "  ";
      }
      Xyce::dout() << std::endl;
    }
  }

  // Because we have saved the inductances from parameters in a local storage,
  //  we need to re-read them just in case the L values were specified as
  // dependent expressions (e.g. dependent on global params that are stepped)
  std::vector< InductorInstanceData* >::iterator
    currentInductor = instanceData.begin();
  std::vector< InductorInstanceData* >::iterator
    endInductor = instanceData.end();

  i=0;
  while( currentInductor != endInductor )
  {
    (*currentInductor)->L = inductorInductances[i];
    (*currentInductor)->baseL = inductorInductances[i];
    ++i;
    ++currentInductor;
  }
  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS)  && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << "processParams: saved inductances:  " << std::endl;
    for (currentInductor = instanceData.begin(), i=0; currentInductor!= endInductor; ++currentInductor, ++i)
    {
      Xyce::dout() << " Inductor " << i << " inductance=" << (*currentInductor)->L << std::endl;
    }
  }

  // set the temperature related stuff.
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
    ++currentData;
  }

  // now that the inductances have changed we need to update the matrix.
  updateInductanceMatrix();

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::updateIntermediateVars
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL, Parallel Computational Sciences
// Creation Date : 03/21/2005
//-----------------------------------------------------------------------------
bool Instance::updateIntermediateVars ()
{
  bool bsuccess = true;

  // This method is never called anymore

  return bsuccess;
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
    inductanceVals[ i ] = ((*currentInductor)->L);
    ++i;
    ++currentInductor;
  }

  // compute the inductance matrix
  for( i=0; i<numInductors; ++i)
  {
    for( int j=0; j<numInductors; ++j)
    {
      LO[i][j] = sqrt( inductanceVals[i]*inductanceVals[j] );
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
  // nothing happens in updateIntermediateVars() anymore
  // bsuccess = updateIntermediateVars ();

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << "Instance::updatePrimaryState------------------" << std::endl
         << "\tname = " << getName() << std::endl;
  }

  double * solVec = extData.nextSolVectorRawPtr;

  // Evaluate the derivatives of all (dependent) coupling coefficients w.r.t
  // their variables.  We need these for both new and old DAE.
  int ncoupcoef=couplingCoefficient.size();
  for (int i=0; i<ncoupcoef; ++i)
  {
    if (expPtrs[i])
    {
      double junk;
      expPtrs[i]->evaluate( junk, couplingCoefficientVarDerivs[i]);
    }
  }

  // get the currents in each inductor
  std::vector< InductorInstanceData* >::iterator
  currentInductor = instanceData.begin();
  std::vector< InductorInstanceData* >::iterator
  endInductor = instanceData.end();
  int i = 0;
  while( currentInductor != endInductor )
  {
    if( (getSolverState().dcopFlag) && ((*currentInductor)->ICGiven) )
    {
      Xyce::dout() << "Applying IC value " << i << " " << (*currentInductor)->IC << std::endl;
      inductorCurrents[ i ] = (*currentInductor)->IC;
    }
    else
    {
      inductorCurrents[ i ] = solVec[ (*currentInductor)->li_Branch ];
    }
    ++i;
    ++currentInductor;
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
// Creator       : Richard Schiek, SNL, Electrical and Microsystems Modeling
// Creation Date : 02/06/07
//-----------------------------------------------------------------------------
void Instance::loadErrorWeightMask ()
{
#ifndef Xyce_NO_MUTINDLIN_MASK
  Linear::Vector * maskVectorPtr = extData.deviceErrorWeightMask_;

  std::vector< InductorInstanceData* >::iterator currentInductor = instanceData.begin();
  std::vector< InductorInstanceData* >::iterator endInductor = instanceData.end();
  while( currentInductor != endInductor )
  {
    (*maskVectorPtr)[((*currentInductor)->li_Branch)] = 0.0;
    ++currentInductor;
  }
#endif
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEQVector
//
// Purpose       : Loads the Q-vector contributions for a single
//                 Mutual Inductor instance.
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
  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << "Instance::loadDAEQVector---------------------------" << std::endl
         << "\tname = " << getName() << std::endl;
  }
  double * qVec = extData.daeQVectorRawPtr;

  // calculate the following product
  // I = column vector of currents
  // L = row vector of inductances
  // LO = matrix = mutualCup * sqrt( L' * L )
  // LOI = column vector = mutualCup * sqrt( L' * L ) * I
  // LOI[1] = mutualCup * sqrt(L[1]*L[1])*I[1]) +
  //          mutualCup * sqrt(L[1]*L[2])*I[2]) + ...
  //          mutualCup * sqrt(L[1]*L[n])*I[n])

  for( int i = 0; i < numInductors; ++i )
  {
    LOI[ i ] = 0.0;
    for( int j = 0; j < numInductors; ++j )
    {
      LOI[i] += mutualCouplingCoef[i][j] * LO[i][j] * inductorCurrents[j];
    }
  }

  // loop over each inductor and load it's Q vector components
  std::vector< InductorInstanceData* >::iterator currentInductor = instanceData.begin();
  std::vector< InductorInstanceData* >::iterator endInductor = instanceData.end();
  int i=0;
  while( currentInductor != endInductor )
  {
    qVec[(*currentInductor)->li_Branch] += LOI[ i ];
    ++i;
    ++currentInductor;
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEFVector
//
// Purpose       : Loads the F-vector contributions for a single
//                 Mutual Inductor instance.
//
// Special Notes :
//
// Scope         : public
// Creator       : Rich Schiek, SNL, Parallel Computational Sciences
// Creation Date : 03/21/2005
//-----------------------------------------------------------------------------
bool Instance::loadDAEFVector ()
{
  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << "Instance::loadDAEFVector---------------------------" << std::endl
         << "\tname = " << getName() << std::endl;
  }

  double * fVec = extData.daeFVectorRawPtr;
  double * solVec = extData.nextSolVectorRawPtr;

  // loop over each inductor and load it's F vector components
  std::vector< InductorInstanceData* >::iterator currentInductor = instanceData.begin();
  std::vector< InductorInstanceData* >::iterator endInductor = instanceData.end();
  int i = 0;

  while( currentInductor != endInductor )
  {
    double current   = solVec[(*currentInductor)->li_Branch];
    double vNodePos  = solVec[(*currentInductor)->li_Pos];
    double vNodeNeg  = solVec[(*currentInductor)->li_Neg];
    fVec[((*currentInductor)->li_Pos)]    +=  scalingRHS * current;
    fVec[((*currentInductor)->li_Neg)]    += -scalingRHS * current;
    fVec[((*currentInductor)->li_Branch)] += -(vNodePos - vNodeNeg);

    if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
    {
      Xyce::dout() << "  Inductor = " << (*currentInductor)->name
           << "\tcurrent = " << current
           << "\tvNodePos = " << vNodePos
           << "\tvNodeNeg = " << vNodeNeg
           << std::endl;
    }

    ++currentInductor;
    ++i;
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEdQdx
//
// Purpose       : Loads the Q-vector contributions for a single
//                 Mutual Inductor instance.
// Scope         : public
// Creator       : Rich Schiek, SNL, Parallel Computational Sciences
// Creation Date : 03/21/2005
//-----------------------------------------------------------------------------
bool Instance::loadDAEdQdx ()
{
  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << "Instance::loadDAEdQdx--------------------------" << std::endl
         << "\tname = " << getName() << std::endl;
  }

  // During the DCOP, the dQdx*pdt matrix is summed into the Jacobian,
  // even though the Q*pdt vector is not summed into the residual.
  //if (!getSolverState().dcopFlag)
  {
    Linear::Matrix & dQdx = *(extData.dQdxMatrixPtr);

    // loop over each inductor and load it's Q vector components
    std::vector< InductorInstanceData* >::iterator currentInductor = instanceData.begin();
    std::vector< InductorInstanceData* >::iterator endInductor = instanceData.end();
    int i=0;
    while( currentInductor != endInductor )
    {
      for( int j=0; j<numInductors; ++j )
      {
        dQdx[((*currentInductor)->li_Branch)]
                     [(*currentInductor)->inductorCurrentOffsets[j]] += mutualCouplingCoef[i][j] * LO[i][j];
      }
      // finally do all the dependent variable terms
      int numdepterms=(*currentInductor)->depVarPairs.size();
      for (int j=0 ; j<numdepterms; ++j)
      {
        int coefficientNumber=(*currentInductor)->depVarPairs[j].first;
        int depVarNumber=(*currentInductor)->depVarPairs[j].second;
        int otherInductor;

        // indexPairs[coefficient] gives the two inductors coupled by that
        // coefficient.  One of them is currentInductor because we saved that
        // coefficient number in our depVarPairs, so the other is the one
        // we need to know

        if (i==indexPairs[coefficientNumber].first)
          otherInductor=indexPairs[coefficientNumber].second;
        else
          otherInductor=indexPairs[coefficientNumber].first;

        dQdx[((*currentInductor)->li_Branch)][(*currentInductor)->ABraEquDepVarOffsets[j]] +=
          couplingCoefficientVarDerivs[coefficientNumber][depVarNumber]
          *LO[i][otherInductor]*inductorCurrents[otherInductor];
      }

      ++i;
      ++currentInductor;
    }
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEdFdx ()
//
// Purpose       : Loads the F-vector contributions for a single
//                 Mutual Inductor instance.
//
// Special Notes :
//
// Scope         : public
// Creator       : Rich Schiek, SNL, Parallel Computational Sciences
// Creation Date : 03/21/2005
//-----------------------------------------------------------------------------
bool Instance::loadDAEdFdx ()
{
  Linear::Matrix & dFdx = *(extData.dFdxMatrixPtr);

  // loop over each inductor and load it's dFdx components
  std::vector< InductorInstanceData* >::iterator currentInductor = instanceData.begin();
  std::vector< InductorInstanceData* >::iterator endInductor = instanceData.end();
  while( currentInductor != endInductor )
  {
    dFdx[((*currentInductor)->li_Pos)]   [((*currentInductor)->APosEquBraVarOffset)]  +=  scalingRHS;
    dFdx[((*currentInductor)->li_Neg)]   [((*currentInductor)->ANegEquBraVarOffset)]  += -scalingRHS;
    dFdx[((*currentInductor)->li_Branch)][((*currentInductor)->ABraEquPosNodeOffset)] += -1.0;
    dFdx[((*currentInductor)->li_Branch)][((*currentInductor)->ABraEquNegNodeOffset)] +=  1.0;

    ++currentInductor;
  }

  return true;
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

  bool bsuccess = true;
  return bsuccess;
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
// Function      : Model::processParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL, Parallel Computational Sciences
// Creation Date : 03/21/2005
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
    tempCoeff1(0.0),
    tempCoeff2(0.0),
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
  os << "Number of MutIndLin instances: " << isize << std::endl;
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


// MutIndLin Master functions:

//-----------------------------------------------------------------------------
// Function      : Master::updateState
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 12/12/08
//-----------------------------------------------------------------------------
bool Master::updateState (double * solVec, double * staVec, double * stoVec)
{
  for (InstanceVector::const_iterator it = getInstanceBegin(); it != getInstanceEnd(); ++it)
  {
     Instance & inst = *(*it);

    // Evaluate the derivatives of all (dependent) coupling coefficients w.r.t
    // their variables.  We need these for both new and old DAE.
    int ncoupcoef=inst.couplingCoefficient.size();
    for (int i=0; i<ncoupcoef; ++i)
    {
      if (inst.expPtrs[i])
      {
        double junk;
        inst.expPtrs[i]->evaluate( junk, inst.couplingCoefficientVarDerivs[i]);
      }
    }

    // get the currents in each inductor
    std::vector< InductorInstanceData* >::iterator
    currentInductor = inst.instanceData.begin();
    std::vector< InductorInstanceData* >::iterator
    endInductor = inst.instanceData.end();
    {
    int i = 0;
    while( currentInductor != endInductor )
    {
      if( (getSolverState().dcopFlag) && ((*currentInductor)->ICGiven) )
      {
        inst.inductorCurrents[ i ] = (*currentInductor)->IC;
      }
      else
      {
        inst.inductorCurrents[ i ] = solVec[ (*currentInductor)->li_Branch ];
      }
      ++i;
      ++currentInductor;
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
// Creation Date : 12/12/08
//-----------------------------------------------------------------------------
bool Master::loadDAEVectors (double * solVec, double * fVec, double *qVec,  double * bVec, double * leadF, double * leadQ, double * junctionV)
{
  for (InstanceVector::const_iterator it = getInstanceBegin(); it != getInstanceEnd(); ++it)
  {
    Instance & inst = *(*it);

    // F-vector:
    // loop over each inductor and load it's F vector components
    std::vector< InductorInstanceData* >::iterator currentInductor = inst.instanceData.begin();
    std::vector< InductorInstanceData* >::iterator endInductor = inst.instanceData.end();

    while( currentInductor != endInductor )
    {
      double current   = solVec[(*currentInductor)->li_Branch];
      double branchCoef = 1.0;
      if( (getSolverState().dcopFlag) &&  ((*currentInductor)->ICGiven))
      {
        current = (*currentInductor)->IC;
        branchCoef = 0.0;
        solVec[(*currentInductor)->li_Branch] = current;
      }
      double vNodePos  = solVec[(*currentInductor)->li_Pos];
      double vNodeNeg  = solVec[(*currentInductor)->li_Neg];

      fVec[((*currentInductor)->li_Pos)]    +=  inst.scalingRHS * current;
      fVec[((*currentInductor)->li_Neg)]    += -inst.scalingRHS * current;
      fVec[((*currentInductor)->li_Branch)] += -(vNodePos - vNodeNeg)*branchCoef;

      if (inst.loadLeadCurrent)
      {
        leadF[(*currentInductor)->li_branch_data] =  inst.scalingRHS * current;       
        junctionV[(*currentInductor)->li_branch_data] = (vNodePos - vNodeNeg);
      }

      ++currentInductor;
    }

    // Q-vector:
    // calculate the following product
    // I = column vector of currents
    // L = row vector of inductances
    // LO = matrix = mutualCup * sqrt( L' * L )
    // LOI = column vector = mutualCup * sqrt( L' * L ) * I
    // LOI[1] = mutualCup * sqrt(L[1]*L[1])*I[1]) +
    //          mutualCup * sqrt(L[1]*L[2])*I[2]) + ...
    //          mutualCup * sqrt(L[1]*L[n])*I[n])

    for( int i = 0; i < inst.numInductors; ++i )
    {
      inst.LOI[ i ] = 0.0;
      for( int j = 0; j < inst.numInductors; ++j )
      {
        inst.LOI[i] += inst.mutualCouplingCoef[i][j] * inst.LO[i][j] * inst.inductorCurrents[j];
      }
    }

    // loop over each inductor and load it's Q vector components
    currentInductor = inst.instanceData.begin();
    endInductor = inst.instanceData.end();
    int li=0;
    while( currentInductor != endInductor )
    {
      qVec[(*currentInductor)->li_Branch] += inst.LOI[ li ];
      ++li;
      ++currentInductor;
    }
  }

  return true;
}

#ifndef Xyce_NONPOINTER_MATRIX_LOAD
//-----------------------------------------------------------------------------
// Function      : Master::loadDAEMatrices
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 12/12/08
//-----------------------------------------------------------------------------
bool Master::loadDAEMatrices (Linear::Matrix & dFdx, Linear::Matrix & dQdx)
{
  for (InstanceVector::const_iterator it = getInstanceBegin(); it != getInstanceEnd(); ++it)
  {
    Instance & inst = *(*it);
    // loop over each inductor and load it's dFdx components
    std::vector< InductorInstanceData* >::iterator currentInductor = inst.instanceData.begin();
    std::vector< InductorInstanceData* >::iterator endInductor = inst.instanceData.end();
    while( currentInductor != endInductor )
    {
    
      if ( getSolverState().dcopFlag && (*currentInductor)->ICGiven )
      {
        // In the case that an initial condition is specified for an
        // inductor, the DC op should be set up like a current source just
        // for the operating point calculation.
        *((*currentInductor)->f_PosEquBraVarPtr)  += 0.0;
        *((*currentInductor)->f_NegEquBraVarPtr)  += 0.0;
        *((*currentInductor)->f_BraEquPosNodePtr) += 0.0;
        *((*currentInductor)->f_BraEquNegNodePtr) += 0.0;
        *((*currentInductor)->f_BraEquBraVarPtr)  += 1.0;
      }
      else
      {
        *((*currentInductor)->f_PosEquBraVarPtr)  +=  inst.scalingRHS;
        *((*currentInductor)->f_NegEquBraVarPtr)  += -inst.scalingRHS;
        *((*currentInductor)->f_BraEquPosNodePtr) += -1.0;
        *((*currentInductor)->f_BraEquNegNodePtr) +=  1.0;
      }

      ++currentInductor;
    }

    // During the DCOP, the dQdx*pdt matrix is summed into the Jacobian,
    // even though the Q*pdt vector is not summed into the residual.
    //if (!getSolverState().dcopFlag)
    {
      // loop over each inductor and load it's Q vector components
      currentInductor = inst.instanceData.begin();
      endInductor = inst.instanceData.end();
      int li=0;
      while( currentInductor != endInductor )
      {
        for( int j=0; j<inst.numInductors; ++j )
        {
          *((*currentInductor)->q_inductorCurrentPtrs[j]) += inst.mutualCouplingCoef[li][j] * inst.LO[li][j];
        }
        // finally do all the dependent variable terms
        int numdepterms=(*currentInductor)->depVarPairs.size();
        for (int j=0 ; j<numdepterms; ++j)
        {
          int coefficientNumber=(*currentInductor)->depVarPairs[j].first;
          int depVarNumber=(*currentInductor)->depVarPairs[j].second;
          int otherInductor;

          // indexPairs[coefficient] gives the two inductors coupled by that
          // coefficient.  One of them is currentInductor because we saved that
          // coefficient number in our depVarPairs, so the other is the one
          // we need to know

          if (li==inst.indexPairs[coefficientNumber].first)
            otherInductor=inst.indexPairs[coefficientNumber].second;
          else
            otherInductor=inst.indexPairs[coefficientNumber].second;

          *((*currentInductor)->q_BraEquDepVarPtrs[j]) +=
                inst.couplingCoefficientVarDerivs[coefficientNumber][depVarNumber]
                *inst.LO[li][otherInductor]*inst.inductorCurrents[otherInductor];
        }

        ++li;
        ++currentInductor;
      }
    }
  }

  return true;
}

#else
//-----------------------------------------------------------------------------
// Function      : Master::loadDAEMatrices
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 12/12/08
//-----------------------------------------------------------------------------
bool Master::loadDAEMatrices (Linear::Matrix & dFdx, Linear::Matrix & dQdx)
{
  for (InstanceVector::const_iterator it = getInstanceBegin(); it != getInstanceEnd(); ++it)
  {
     Instance & inst = *(*it);

    // loop over each inductor and load it's dFdx components
    std::vector< InductorInstanceData* >::iterator currentInductor = inst.instanceData.begin();
    std::vector< InductorInstanceData* >::iterator endInductor = inst.instanceData.end();
    while( currentInductor != endInductor )
    {
      dFdx[((*currentInductor)->li_Pos)]   [((*currentInductor)->APosEquBraVarOffset)]  +=  inst.scalingRHS;
      dFdx[((*currentInductor)->li_Neg)]   [((*currentInductor)->ANegEquBraVarOffset)]  += -inst.scalingRHS;
      dFdx[((*currentInductor)->li_Branch)][((*currentInductor)->ABraEquPosNodeOffset)] += -1.0;
      dFdx[((*currentInductor)->li_Branch)][((*currentInductor)->ABraEquNegNodeOffset)] +=  1.0;

      ++currentInductor;
    }

    // During the DCOP, the dQdx*pdt matrix is summed into the Jacobian,
    // even though the Q*pdt vector is not summed into the residual.
    //if (!getSolverState().dcopFlag)
    {
      // loop over each inductor and load it's Q vector components
      currentInductor = inst.instanceData.begin();
      endInductor = inst.instanceData.end();
      int i=0;
      while( currentInductor != endInductor )
      {
        for( int j=0; j<inst.numInductors; ++j )
        {
          dQdx[((*currentInductor)->li_Branch)]
                      [(*currentInductor)->inductorCurrentOffsets[j]] += inst.mutualCouplingCoef[i][j] * inst.LO[i][j];
        }
        // finally do all the dependent variable terms
        int numdepterms=(*currentInductor)->depVarPairs.size();
        for (int j=0 ; j<numdepterms; ++j)
        {
          int coefficientNumber=(*currentInductor)->depVarPairs[j].first;
          int depVarNumber=(*currentInductor)->depVarPairs[j].second;
          int otherInductor;

          // indexPairs[coefficient] gives the two inductors coupled by that
          // coefficient.  One of them is currentInductor because we saved that
          // coefficient number in our depVarPairs, so the other is the one
          // we need to know

          if (i==inst.indexPairs[coefficientNumber].first)
            otherInductor=inst.indexPairs[coefficientNumber].second;
          else
            otherInductor=inst.indexPairs[coefficientNumber].second;

          dQdx[((*currentInductor)->li_Branch)][(*currentInductor)->ABraEquDepVarOffsets[j]] +=
            inst.couplingCoefficientVarDerivs[coefficientNumber][depVarNumber]
            *inst.LO[i][otherInductor]*inst.inductorCurrents[otherInductor];
        }

        ++i;
        ++currentInductor;
      }
    }
  }
  return true;
}
#endif

Device *Traits::factory(const Configuration &configuration, const FactoryBlock &factory_block)
{

  return new Master(configuration, factory_block, factory_block.solverState_, factory_block.deviceOptions_);
}

void registerDevice()
{
  Config<Traits>::addConfiguration()
    .registerDevice("mil", 1)
    .registerModelType("mil", 1);
}

} // namespace MutIndLin
} // namespace Device
} // namespace Xyce

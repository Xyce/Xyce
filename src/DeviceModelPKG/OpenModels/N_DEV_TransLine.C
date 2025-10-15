//-------------------------------------------------------------------------
//   Copyright 2002-2025 National Technology & Engineering Solutions of
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
// Creator        : Eric Keiter, SNL
//
// Creation Date  : 9/17/2013
//
//
//
//
//-------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------
#include <algorithm>
#include <N_UTL_Math.h>

// ----------   Xyce Includes   ----------
#include <N_DEV_Const.h>
#include <N_DEV_DeviceMaster.h>
#include <N_DEV_DeviceOptions.h>
#include <N_DEV_ExternData.h>
#include <N_DEV_MatrixLoadData.h>
#include <N_DEV_Message.h>
#include <N_DEV_SolverState.h>
#include <N_DEV_TransLine.h>
#include <N_ERH_ErrorMgr.h>
#include <N_LAS_Matrix.h>
#include <N_LAS_Vector.h>
#include <N_UTL_Expression.h>
#include <N_UTL_FeatureTest.h>

#include <Teuchos_BLAS.hpp>
#include <Teuchos_Utils.hpp>
#include <Teuchos_LAPACK.hpp>

#undef HAVE_INTTYPES_H
#undef HAVE_STDINT_H
#include <Trilinos_Util.h>

namespace Xyce {
namespace Device {

namespace TransLine {

void Traits::loadInstanceParameters(ParametricData<TransLine::Instance> &p)
{
  p.addPar ("LUMPS",1,&TransLine::Instance::numLumps)
   .setGivenMember(&TransLine::Instance::numLumpsGiven);

  p.addPar ("LEN",0.0,&TransLine::Instance::length)
   .setGivenMember(&TransLine::Instance::lengthGiven)
   .setUnit(U_METER)
   .setCategory(CAT_NONE)
   .setDescription("length of line");
}

void Traits::loadModelParameters(ParametricData<TransLine::Model> &p)
{
  p.addPar ("R",0.0,&TransLine::Model::resist)
   .setGivenMember(&TransLine::Model::resistGiven)
   .setUnit(U_OHMMM1)
   .setCategory(CAT_NONE)
   .setDescription("Resistance per unit length");

  p.addPar ("L",0.0,&TransLine::Model::induct)
   .setGivenMember(&TransLine::Model::inductGiven)
   .setUnit(U_HMM1)
   .setCategory(CAT_NONE)
   .setDescription("Inductance per unit length");

  p.addPar ("G",0.0,&TransLine::Model::conduct)
   .setGivenMember(&TransLine::Model::conductGiven)
   .setUnit(U_OHMM1MM1)
   .setCategory(CAT_NONE)
   .setDescription("Conductance per unit length");

  p.addPar ("C",0.0,&TransLine::Model::capac)
   .setGivenMember(&TransLine::Model::capacGiven)
   .setUnit(U_FARADMM1)
   .setCategory(CAT_NONE)
   .setDescription("Capacitance per unit length");

  p.addPar ("ELEV",2,&TransLine::Model::elevNumber)
    .setGivenMember(&TransLine::Model::elevNumberGiven);
}

// Class Instance
//-----------------------------------------------------------------------------
// Function      : Instance::processParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 9/17/2013
//-----------------------------------------------------------------------------
bool Instance::processParams ()
{
  // setup R,L,C,G.
  double dx = (1.0/static_cast<double>(numLumps))*length;
  if (model_.resistGiven && (model_.resist != 0.0))
  {
    R = model_.resist*dx;
    G = 1/R;
  }
  if (model_.capacGiven)
  {
    C = model_.capac*dx;
  }
  if (model_.inductGiven)
  {
    L = model_.induct*dx;
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::Instance
// Purpose       : instance block constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 9/17/2013
//-----------------------------------------------------------------------------
Instance::Instance(
  const Configuration & configuration,
  const InstanceBlock & IB,
  Model & TransLineiter,
  const FactoryBlock &  factory_block)
  : DeviceInstance(IB, configuration.getInstanceParameters(), factory_block),
    model_(TransLineiter),
    numLumps(1),
    numTransLineVars(0),
    length(0.0),
    numLumpsGiven(false),
    lengthGiven(false),
    L(0.0),C(0.0),G(0.0),R(0.0)
{
  numExtVars   = IB.numExtVars;  // we have as many as were specified on the
                                 // instance line

  // Set params to constant default values:
  setDefaultParams ();

  // Set params according to instance line and constant defaults from metadata:
  setParams (IB.params);

  // Calculate any parameters specified as expressions:
  updateDependentParameters();

  // calculate dependent (ie computed) params:
  processParams ();

  // Initialize which case this is based on the nonzero user-specified
  // parameters.
  if ((model_.resist == 0) && (model_.conduct == 0) &&
      (model_.capac != 0) && (model_.induct != 0))
  {
    model_.specialCase = TRANS_MOD_LC;
  }
  else if ((model_.resist != 0) && (model_.conduct == 0) &&
           (model_.capac != 0) && (model_.induct != 0))
  {
    model_.specialCase = TRANS_MOD_RLC;
  }

  if (DEBUG_DEVICE)
  {
    if (isActive(Diag::DEVICE_PARAMETERS))
    {
      std::cout << "model_.resist = " << model_.resist <<std::endl;
      std::cout << "model_.capac  = " << model_.capac <<std::endl;
      std::cout << "model_.induct = " << model_.induct <<std::endl;

      std::cout << "R = " << R <<std::endl;
      std::cout << "G = " << G <<std::endl;
      std::cout << "C = " << C <<std::endl;
      std::cout << "L = " << L <<std::endl;

      if (model_.specialCase == TRANS_MOD_RLC)
      {
        std::cout << "RLC line" <<std::endl;
      }
      else if (model_.specialCase == TRANS_MOD_LC)
      {
        std::cout << "LC line" <<std::endl;
      }
    }
  }

  if ( model_.specialCase == TRANS_MOD_RLC)
  {
    // RLC line.
    //
    // each RLC lump has 2 voltage nodes and 1 branch current (3 vars)
    // there are two external variables, at opposite ends of the line.
    // One of them has to be subtracted from the numIntVars.
    //
    //           lump 1        lump 2       lump 3        lump 4
    //  Input --L--1--R--2  --L--3--R--4  --L--5--R--6  --L--7--R--x Output
    //                  |             |             |             |
    //                  C             C             C             C
    //                  |             |             |             |
    //                  0(gnd)        0(gnd)        0(gnd)        0(gnd)
    //
    //
    numIntVars   = numLumps*3-1;
    numExtVars   = 2;
    numStateVars = 0;

    lumpVec.resize(numLumps);

    // loop over the lumps, set up nominal indices
    for (int lump=0;lump<numLumps;++lump)
    {
      int varsPerLump=3;  // 2 voltage nodes and one branch current.
      lumpVec[lump].indexV1 = lump*varsPerLump;   // voltage node 1
      lumpVec[lump].indexV2 = lump*varsPerLump+1; // voltage node 2
      lumpVec[lump].indexI  = lump*varsPerLump+2; // branch current
      lumpVec[lump].indexV3 = (lump+1)*varsPerLump; // voltage node 3
    }

    if (DEBUG_DEVICE)
    {
      if (isActive(Diag::DEVICE_PARAMETERS))
      {
        std::cout << "lumps = " << numLumps <<std::endl;
        for (int lump=0;lump<numLumps;++lump)
        {
          std::cout << "lumpVec["<<lump<<"]: v1 = " << lumpVec[lump].indexV1;
          std::cout << "  v2 = " << lumpVec[lump].indexV2;
          std::cout << "  i  = " << lumpVec[lump].indexI ;
          std::cout << "  v3 = " << lumpVec[lump].indexV3;
          std::cout << std::endl;
        }
      }
    }

    // now correct the lump vector to account for the ends of the
    // line, which are external vars.  The first node (input) will not cause a
    // correction as it is naturally the first index(0) anyway.
    //
    // All nodes that are greater than 0 will need to be incremented by 1,
    // as the next external node will be index 1.
    for (int lump=0;lump<numLumps;++lump)
    {
      if (lumpVec[lump].indexV1 > 0) lumpVec[lump].indexV1++;
      if (lumpVec[lump].indexV2 > 0) lumpVec[lump].indexV2++;
      if (lumpVec[lump].indexI  > 0) lumpVec[lump].indexI++;
      if (lumpVec[lump].indexV3 > 0) lumpVec[lump].indexV3++;
    }

    // the last index (being an external var) much be set back to 1.
    lumpVec[numLumps-1].indexV3 = 1;


    if (DEBUG_DEVICE)
    {
      if (isActive(Diag::DEVICE_PARAMETERS))
      {
        std::cout << "lumps = " << numLumps <<std::endl;
        for (int lump=0;lump<numLumps;++lump)
        {
          std::cout << "lumpVec["<<lump<<"]: v1 = " << lumpVec[lump].indexV1;
          std::cout << "  v2 = " << lumpVec[lump].indexV2;
          std::cout << "  i  = " << lumpVec[lump].indexI ;
          std::cout << "  v3 = " << lumpVec[lump].indexV3;
          std::cout << std::endl;
        }
      }
    }

    // setup jacobian stamp
    jacStamp.resize(numExtVars+numIntVars);

    // do jacStamp for the exterior lumps
    {
      int lump=0;
      int n1 = lumpVec[lump].indexV1;
      int n2 = lumpVec[lump].indexV2;
      int ii = lumpVec[lump].indexI;
      int n3 = lumpVec[lump].indexV3;

      jacStamp[n1].resize(1);
      jacStamp[n1][0] = ii;
      jacStamp[n2].resize(3);
      jacStamp[n2][0] = n2;
      jacStamp[n2][1] = ii;
      jacStamp[n2][2] = n3;
      jacStamp[ii].resize(3);
      jacStamp[ii][0] = n1;
      jacStamp[ii][1] = n2;
      jacStamp[ii][2] = ii;

      if (numLumps==1)
      {
        jacStamp[n3].resize(2);
        jacStamp[n3][0] = n2;
        jacStamp[n3][1] = n3;
      }
    }

    if (numLumps>1)
    {
      int lump=numLumps-1;
      int n1 = lumpVec[lump].indexV1;
      int n2 = lumpVec[lump].indexV2;
      int ii = lumpVec[lump].indexI;
      int n3 = lumpVec[lump].indexV3;

      jacStamp[n1].resize(3);
      jacStamp[n1][0] = lumpVec[lump-1].indexV2;
      jacStamp[n1][1] = n1;
      jacStamp[n1][2] = ii;

      jacStamp[n2].resize(3);
      jacStamp[n2][0] = n2;
      jacStamp[n2][1] = ii;
      jacStamp[n2][2] = n3;

      jacStamp[ii].resize(3);
      jacStamp[ii][0] = n1;
      jacStamp[ii][1] = n2;
      jacStamp[ii][2] = ii;

      jacStamp[n3].resize(2);
      jacStamp[n3][0] = n2;
      jacStamp[n3][1] = n3;

    }

    // do jacStamp for the interior lumps
    for (int lump=1;lump<numLumps-1;++lump)
    {
      int n1 = lumpVec[lump].indexV1;
      int n2 = lumpVec[lump].indexV2;
      int ii = lumpVec[lump].indexI;
      int n3 = lumpVec[lump].indexV3;

      jacStamp[n1].resize(3);
      jacStamp[n1][0] = lumpVec[lump-1].indexV2;
      jacStamp[n1][1] = n1;
      jacStamp[n1][2] = ii;
      jacStamp[n2].resize(3);
      jacStamp[n2][0] = n2;
      jacStamp[n2][1] = ii;
      jacStamp[n2][2] = n3;
      jacStamp[ii].resize(3);
      jacStamp[ii][0] = n1;
      jacStamp[ii][1] = n2;
      jacStamp[ii][2] = ii;
    }
  }
  else if ( model_.specialCase == TRANS_MOD_LC)
  {
    // LC line.
    //
    // each LC lump has 1 voltage node and 1 branch current (2 vars)
    // there are two external variables, at opposite ends of the line.
    // One of them has to be subtracted from the numIntVars.
    //
    //           lump 1        lump 2       lump 3        lump 4
    //  Input --L--------1----L--------2----L--------3----L--------x Output
    //                  |             |             |             |
    //                  C             C             C             C
    //                  |             |             |             |
    //                  0(gnd)        0(gnd)        0(gnd)        0(gnd)
    //
    numIntVars   = numLumps*2-1;
    numExtVars   = 2;
    numStateVars = 0;

    lumpVec.resize(numLumps);

    // loop over the lumps, set up nominal indices
    for (int lump=0;lump<numLumps;++lump)
    {
      int varsPerLump=2;  // 2 voltage nodes and one branch current.
      lumpVec[lump].indexV1 = lump*varsPerLump;   // voltage node 1
      lumpVec[lump].indexI  = lump*varsPerLump+1; // branch current
      lumpVec[lump].indexV2 = (lump+1)*varsPerLump; // voltage node 3
    }

    if (DEBUG_DEVICE)
    {
      if (isActive(Diag::DEVICE_PARAMETERS))
      {
        std::cout << "lumps = " << numLumps <<std::endl;
        for (int lump=0;lump<numLumps;++lump)
        {
          std::cout << "lumpVec["<<lump<<"]: v1 = " << lumpVec[lump].indexV1;
          std::cout << "  i  = " << lumpVec[lump].indexI ;
          std::cout << "  v2 = " << lumpVec[lump].indexV2;
          std::cout << std::endl;
        }
      }
    }

    // now correct the lump vector to account for the ends of the
    // line, which are external vars.  The first node (input) will not cause a
    // correction as it is naturally the first index(0) anyway.
    //
    // All nodes that are greater than 0 will need to be incremented by 1,
    // as the next external node will be index 1.
    for (int lump=0;lump<numLumps;++lump)
    {
      if (lumpVec[lump].indexV1 > 0) lumpVec[lump].indexV1++;
      if (lumpVec[lump].indexI  > 0) lumpVec[lump].indexI++;
      if (lumpVec[lump].indexV2 > 0) lumpVec[lump].indexV2++;
    }

    // the last index (being an external var) much be set back to 1.
    lumpVec[numLumps-1].indexV2 = 1;


    if (DEBUG_DEVICE)
    {
      if (isActive(Diag::DEVICE_PARAMETERS))
      {
        std::cout << "lumps = " << numLumps <<std::endl;
        for (int lump=0;lump<numLumps;++lump)
        {
          std::cout << "lumpVec["<<lump<<"]: v1 = " << lumpVec[lump].indexV1;
          std::cout << "  i  = " << lumpVec[lump].indexI ;
          std::cout << "  v2 = " << lumpVec[lump].indexV2;
          std::cout << std::endl;
        }
      }
    }

    // setup jacobian stamp
    jacStamp.resize(numExtVars+numIntVars);

    // do jacStamp for the exterior lumps
    {
      int lump=0;
      int n1 = lumpVec[lump].indexV1;
      int n2 = lumpVec[lump].indexV2;
      int ii = lumpVec[lump].indexI;

      jacStamp[n1].resize(1);
      jacStamp[n1][0] = ii;

      jacStamp[ii].resize(3);
      jacStamp[ii][0] = n1;
      jacStamp[ii][1] = ii;
      jacStamp[ii][2] = n2;

      if (numLumps==1)
      {
        jacStamp[n2].resize(2);
        jacStamp[n2][0] = ii;
        jacStamp[n2][1] = n2;
      }
    }

    if (numLumps>1)
    {
      int lump=numLumps-1;
      int n1 = lumpVec[lump].indexV1;
      int n2 = lumpVec[lump].indexV2;
      int ii = lumpVec[lump].indexI;

      jacStamp[n1].resize(3);
      jacStamp[n1][0] = lumpVec[lump-1].indexI;
      jacStamp[n1][1] = n1;
      jacStamp[n1][2] = ii;

      jacStamp[ii].resize(3);
      jacStamp[ii][0] = n1;
      jacStamp[ii][1] = ii;
      jacStamp[ii][2] = n2;

      jacStamp[n2].resize(2);
      jacStamp[n2][0] = ii;
      jacStamp[n2][1] = n2;
    }

    // do jacStamp for the interior lumps
    for (int lump=1;lump<numLumps-1;++lump)
    {
      int n1 = lumpVec[lump].indexV1;
      int n2 = lumpVec[lump].indexV2;
      int ii = lumpVec[lump].indexI;

      jacStamp[n1].resize(3);
      jacStamp[n1][0] = lumpVec[lump-1].indexI;
      jacStamp[n1][1] = n1;
      jacStamp[n1][2] = ii;

      jacStamp[ii].resize(3);
      jacStamp[ii][0] = n1;
      jacStamp[ii][1] = ii;
      jacStamp[ii][2] = n2;
    }
  }

  if (DEBUG_DEVICE)
  {
    if (isActive(Diag::DEVICE_PARAMETERS))
    {
      int size=jacStamp.size();
      for (int i=0; i<size; ++i)
      {
        int sizeJ = jacStamp[i].size();
        for (int j=0;j<sizeJ;++j)
        {
          std::cout << "jacStamp["<<i<<"]["<<j<<"] = " << jacStamp[i][j];
          std::cout << std::endl;
        }
      }
    }

  }

}

//-----------------------------------------------------------------------------
// Function      : Instance::~Instance
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 9/17/2013
//-----------------------------------------------------------------------------
Instance::~Instance()
{
}

// Additional Declarations

//-----------------------------------------------------------------------------
// Function      : Instance::registerLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 9/17/2013
//-----------------------------------------------------------------------------
void Instance::registerLIDs( const std::vector<int> & intLIDVecRef,
                             const std::vector<int> & extLIDVecRef)
{
  AssertLIDs(intLIDVecRef.size() == numIntVars);
  AssertLIDs(extLIDVecRef.size() == numExtVars);

  // Copy over the local ID lists:
  intLIDVec = intLIDVecRef;
  extLIDVec = extLIDVecRef;

  // Now use these lists to obtain the indices into the linear algebra entities.

  if ( model_.specialCase == TRANS_MOD_RLC)
  {
    lumpVec[0].li_V1 = extLIDVec[0];
    lumpVec[numLumps-1].li_V3 = extLIDVec[1];

    int lid=0;
    lumpVec[0].li_V2 = intLIDVec[lid++];
    lumpVec[0].li_I  = intLIDVec[lid++];

    if (numLumps>1)
    {
      for (int i=1;i<numLumps-1;++i)
      {
        lumpVec[i].li_V1 =  intLIDVec[lid++];
        lumpVec[i].li_V2 =  intLIDVec[lid++];
        lumpVec[i].li_I  =  intLIDVec[lid++];
      }

      lumpVec[numLumps-1].li_V1 = intLIDVec[lid++];
      lumpVec[numLumps-1].li_V2 = intLIDVec[lid++];
      lumpVec[numLumps-1].li_I  = intLIDVec[lid++];

      for (int i=0;i<numLumps-1;++i)
      {
        lumpVec[i].li_V3 = lumpVec[i+1].li_V1;
      }
    }
  }
  else if ( model_.specialCase == TRANS_MOD_LC)
  {
    lumpVec[0].li_V1 = extLIDVec[0];
    lumpVec[numLumps-1].li_V2 = extLIDVec[1];

    int lid=0;
    lumpVec[0].li_I  = intLIDVec[lid++];

    if (numLumps>1)
    {
      for (int i=1;i<numLumps-1;++i)
      {
        lumpVec[i].li_V1 =  intLIDVec[lid++];
        lumpVec[i].li_I  =  intLIDVec[lid++];
      }

      lumpVec[numLumps-1].li_V1 = intLIDVec[lid++];
      lumpVec[numLumps-1].li_I  = intLIDVec[lid++];

      for (int i=0;i<numLumps-1;++i)
      {
        lumpVec[i].li_V2 = lumpVec[i+1].li_V1;
      }
    }
  }

  if (DEBUG_DEVICE)
  {
    if (isActive(Diag::DEVICE_PARAMETERS))
    {
      for (int i=0; i<numIntVars; i++)
      {
        std::cout << "intLIDVec["<<i<<"] = " << intLIDVec[i] << std::endl;
      }

      if ( model_.specialCase == TRANS_MOD_RLC)
      {
        for (int i=0; i<numLumps; i++)
        {
          std::cout << "lumpVec["<<i<<"].li_V1 = " << lumpVec[i].li_V1 <<std::endl;
          std::cout << "lumpVec["<<i<<"].li_V2 = " << lumpVec[i].li_V2 <<std::endl;
          std::cout << "lumpVec["<<i<<"].li_I  = " << lumpVec[i].li_I  <<std::endl;
          std::cout << "lumpVec["<<i<<"].li_V3 = " << lumpVec[i].li_V3 <<std::endl;
        }
      }
      else if ( model_.specialCase == TRANS_MOD_LC)
      {
        for (int i=0; i<numLumps; i++)
        {
          std::cout << "lumpVec["<<i<<"].li_V1 = " << lumpVec[i].li_V1 <<std::endl;
          std::cout << "lumpVec["<<i<<"].li_I  = " << lumpVec[i].li_I  <<std::endl;
          std::cout << "lumpVec["<<i<<"].li_V2 = " << lumpVec[i].li_V2 <<std::endl;
        }
      }
    }
  }

}

//-----------------------------------------------------------------------------
// Function      : Instance::registerStateLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 9/17/2013
//-----------------------------------------------------------------------------
void Instance::registerStateLIDs( const std::vector<int> & staLIDVecRef)
{
  AssertLIDs(staLIDVecRef.size() == numStateVars);

  // Copy over the global ID lists:
  staLIDVec = staLIDVecRef;
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
  for (int i = 0; i < numLumps; i++)
  {
    if (model_.specialCase == TRANS_MOD_RLC)
    {
      if (i != 0)
        addInternalNode(symbol_table, lumpVec[i].li_V1, getName(), "v1_lump" + Teuchos::Utils::toString(i));
      addInternalNode(symbol_table, lumpVec[i].li_V2, getName(), "v2_lump" + Teuchos::Utils::toString(i));
      addInternalNode(symbol_table, lumpVec[i].li_I, getName(), "I__lump" + Teuchos::Utils::toString(i));
    }
    else if (model_.specialCase == TRANS_MOD_LC)
    {
      if (i != 0)
        addInternalNode(symbol_table, lumpVec[i].li_V1, getName(), "v1_lump" + Teuchos::Utils::toString(i));
      addInternalNode(symbol_table, lumpVec[i].li_I, getName(), "I__lump" + Teuchos::Utils::toString(i));
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : Instance::jacobianStamp
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 9/17/2013
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
// Creator       : Eric Keiter, SNL
// Creation Date : 9/17/2013
//-----------------------------------------------------------------------------
void Instance::registerJacLIDs( const std::vector< std::vector<int> > & jacLIDVec )
{
  DeviceInstance::registerJacLIDs( jacLIDVec );


  if ( model_.specialCase == TRANS_MOD_RLC)
  {
    // do jacStamp for the exterior lumps
    {
      int lump=0;
      int n1 = lumpVec[lump].indexV1;
      int n2 = lumpVec[lump].indexV2;
      int ii = lumpVec[lump].indexI;
      int n3 = lumpVec[lump].indexV3;

      lumpVec[lump].offset_v1_ii = jacLIDVec[n1][0];

      lumpVec[lump].offset_v2_v2 = jacLIDVec[n2][0];
      lumpVec[lump].offset_v2_ii = jacLIDVec[n2][1];
      lumpVec[lump].offset_v2_v3 = jacLIDVec[n2][2];

      lumpVec[lump].offset_ii_v1 = jacLIDVec[ii][0];
      lumpVec[lump].offset_ii_v2 = jacLIDVec[ii][1];
      lumpVec[lump].offset_ii_ii = jacLIDVec[ii][2];

      if (numLumps==1)
      {
        lumpVec[lump].offset_v3_v2 = jacLIDVec[n3][0];
        lumpVec[lump].offset_v3_v3 = jacLIDVec[n3][1];
      }
    }

    if (numLumps>1)
    {
      int lump=numLumps-1;
      int n1 = lumpVec[lump].indexV1;
      int n2 = lumpVec[lump].indexV2;
      int ii = lumpVec[lump].indexI;
      int n3 = lumpVec[lump].indexV3;

      lumpVec[lump].offset_v1_v2m1 = jacLIDVec[n1][0];
      lumpVec[lump].offset_v1_v1   = jacLIDVec[n1][1];
      lumpVec[lump].offset_v1_ii   = jacLIDVec[n1][2];

      lumpVec[lump].offset_v2_v2   = jacLIDVec[n2][0];
      lumpVec[lump].offset_v2_ii   = jacLIDVec[n2][1];
      lumpVec[lump].offset_v2_v3   = jacLIDVec[n2][2];

      lumpVec[lump].offset_ii_v1   = jacLIDVec[ii][0];
      lumpVec[lump].offset_ii_v2   = jacLIDVec[ii][1];
      lumpVec[lump].offset_ii_ii   = jacLIDVec[ii][2];

      lumpVec[lump].offset_v3_v2   = jacLIDVec[n3][0];
      lumpVec[lump].offset_v3_v3   = jacLIDVec[n3][1];
    }

    // do jacLIDVec for the interior lumps
    for (int lump=1;lump<numLumps-1;++lump)
    {
      int n1 = lumpVec[lump].indexV1;
      int n2 = lumpVec[lump].indexV2;
      int ii = lumpVec[lump].indexI;
      int n3 = lumpVec[lump].indexV3;

      lumpVec[lump].offset_v1_v2m1 = jacLIDVec[n1][0];
      lumpVec[lump].offset_v1_v1   = jacLIDVec[n1][1];
      lumpVec[lump].offset_v1_ii   = jacLIDVec[n1][2];

      lumpVec[lump].offset_v2_v2   = jacLIDVec[n2][0];
      lumpVec[lump].offset_v2_ii   = jacLIDVec[n2][1];
      lumpVec[lump].offset_v2_v3   = jacLIDVec[n2][2];

      lumpVec[lump].offset_ii_v1   = jacLIDVec[ii][0];
      lumpVec[lump].offset_ii_v2   = jacLIDVec[ii][1];
      lumpVec[lump].offset_ii_ii   = jacLIDVec[ii][2];
    }

    // add the v3 offsets
    for (int lump=0; lump<numLumps-1; ++lump)
    {
      lumpVec[lump].offset_v3_v2   = lumpVec[lump+1].offset_v1_v2m1;
      lumpVec[lump].offset_v3_v3   = lumpVec[lump+1].offset_v1_v1;
    }
  }
  else if ( model_.specialCase == TRANS_MOD_LC)
  {
  // do jacStamp for the exterior lumps
    {
      int lump=0;
      int n1 = lumpVec[lump].indexV1;
      int n2 = lumpVec[lump].indexV2;
      int ii = lumpVec[lump].indexI;

      lumpVec[lump].offset_v1_ii = jacLIDVec[n1][0];

      lumpVec[lump].offset_ii_v1 = jacLIDVec[ii][0];
      lumpVec[lump].offset_ii_ii = jacLIDVec[ii][1];
      lumpVec[lump].offset_ii_v2 = jacLIDVec[ii][2];

      if (numLumps==1)
      {
        lumpVec[lump].offset_v2_ii = jacLIDVec[n2][0];
        lumpVec[lump].offset_v2_v2 = jacLIDVec[n2][1];
      }
    }

    if (numLumps>1)
    {
      int lump=numLumps-1;
      int n1 = lumpVec[lump].indexV1;
      int n2 = lumpVec[lump].indexV2;
      int ii = lumpVec[lump].indexI;

      lumpVec[lump].offset_v1_iim1 = jacLIDVec[n1][0];
      lumpVec[lump].offset_v1_v1   = jacLIDVec[n1][1];
      lumpVec[lump].offset_v1_ii   = jacLIDVec[n1][2];

      lumpVec[lump].offset_ii_v1   = jacLIDVec[ii][0];
      lumpVec[lump].offset_ii_ii   = jacLIDVec[ii][1];
      lumpVec[lump].offset_ii_v2   = jacLIDVec[ii][2];

      lumpVec[lump].offset_v2_ii   = jacLIDVec[n2][0];
      lumpVec[lump].offset_v2_v2   = jacLIDVec[n2][1];
    }

    // do jacLIDVec for the interior lumps
    for (int lump=1;lump<numLumps-1;++lump)
    {
      int n1 = lumpVec[lump].indexV1;
      int n2 = lumpVec[lump].indexV2;
      int ii = lumpVec[lump].indexI;

      lumpVec[lump].offset_v1_iim1 = jacLIDVec[n1][0];
      lumpVec[lump].offset_v1_v1   = jacLIDVec[n1][1];
      lumpVec[lump].offset_v1_ii   = jacLIDVec[n1][2];

      lumpVec[lump].offset_ii_v1   = jacLIDVec[ii][0];
      lumpVec[lump].offset_ii_ii   = jacLIDVec[ii][1];
      lumpVec[lump].offset_ii_v2   = jacLIDVec[ii][2];
    }

    // add the v2 offsets
    for (int lump=0; lump<numLumps-1; ++lump)
    {
      lumpVec[lump].offset_v2_ii   = lumpVec[lump+1].offset_v1_iim1;
      lumpVec[lump].offset_v2_v2   = lumpVec[lump+1].offset_v1_v1;
    }
  }


  if (DEBUG_DEVICE)
  {
    if (isActive(Diag::DEVICE_PARAMETERS))
    {
      if ( model_.specialCase == TRANS_MOD_RLC)
      {
        for (int i=0;i<numLumps;++i)
        {
          std::cout << "lump = " << i <<std::endl;

          std::cout << "offset_v1_v2m1 = " << lumpVec[i].offset_v1_v2m1 << std::endl;
          std::cout << "offset_v1_v1   = " << lumpVec[i].offset_v1_v1 << std::endl;
          std::cout << "offset_v1_ii   = " << lumpVec[i].offset_v1_ii << std::endl;


          std::cout << "offset_v2_v2   = " << lumpVec[i].offset_v2_v2 << std::endl;
          std::cout << "offset_v2_ii   = " << lumpVec[i].offset_v2_ii << std::endl;
          std::cout << "offset_v2_v3   = " << lumpVec[i].offset_v2_v3 << std::endl;


          std::cout << "offset_ii_v1   = " << lumpVec[i].offset_ii_v1 << std::endl;
          std::cout << "offset_ii_v2   = " << lumpVec[i].offset_ii_v2 << std::endl;
          std::cout << "offset_ii_ii   = " << lumpVec[i].offset_ii_ii << std::endl;

          std::cout << "offset_v3_v2   = " << lumpVec[i].offset_v3_v2 << std::endl;
          std::cout << "offset_v3_v3   = " << lumpVec[i].offset_v3_v3 << std::endl;
        }
      }
      else if ( model_.specialCase == TRANS_MOD_LC)
      {
        for (int i=0;i<numLumps;++i)
        {
          std::cout << "lump = " << i <<std::endl;

          std::cout << "offset_v1_iim1 = " << lumpVec[i].offset_v1_iim1 << std::endl;
          std::cout << "offset_v1_v1   = " << lumpVec[i].offset_v1_v1 << std::endl;
          std::cout << "offset_v1_ii   = " << lumpVec[i].offset_v1_ii << std::endl;

          std::cout << "offset_ii_v1   = " << lumpVec[i].offset_ii_v1 << std::endl;
          std::cout << "offset_ii_ii   = " << lumpVec[i].offset_ii_ii << std::endl;
          std::cout << "offset_ii_v2   = " << lumpVec[i].offset_ii_v2 << std::endl;

          std::cout << "offset_v2_ii   = " << lumpVec[i].offset_v2_ii << std::endl;
          std::cout << "offset_v2_v2   = " << lumpVec[i].offset_v2_v2 << std::endl;
       }
      }
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : Instance::setupPointers
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 9/17/2013
//-----------------------------------------------------------------------------
void Instance::setupPointers ()
{

#ifndef Xyce_NONPOINTER_MATRIX_LOAD
#endif
}


//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEQVector
//
// Purpose       : Loads the Q-vector contributions for a single
//                 TransLine instance.
//
// Special Notes : The "Q" std::vector is part of a standard DAE formalism in
//                 which the system of equations is represented as:
//
//                 f(x) = dQ(x)/dt + F(x) - B(t) = 0
//
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 9/17/2013
//-----------------------------------------------------------------------------
bool Instance::loadDAEQVector ()
{
  double * qVec = extData.daeQVectorRawPtr;
  double * solVec = extData.nextSolVectorRawPtr;

  if ( model_.specialCase == TRANS_MOD_RLC)
  {
    for (int i=0; i<numLumps;++i)
    {
      qVec[lumpVec[i].li_I] += L*solVec[lumpVec[i].li_I];
      qVec[lumpVec[i].li_V3] += C*solVec[lumpVec[i].li_V3];
    }
  }
  else if ( model_.specialCase == TRANS_MOD_LC)
  {
    for (int i=0; i<numLumps;++i)
    {
      qVec[lumpVec[i].li_I] += L*solVec[lumpVec[i].li_I];
      qVec[lumpVec[i].li_V2] += C*solVec[lumpVec[i].li_V2];
    }
  }
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEFVector
//
// Purpose       : Loads the F-vector contributions for a single
//                 TransLine instance.
//
// Special Notes :
//
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 9/17/2013
//-----------------------------------------------------------------------------
bool Instance::loadDAEFVector ()
{
  double * fVec = extData.daeFVectorRawPtr;
  double * solVec = extData.nextSolVectorRawPtr;

  if ( model_.specialCase == TRANS_MOD_RLC)
  {
    for (int i=0; i<numLumps; i++)
    {
      // inductor
      fVec[lumpVec[i].li_V1]  += solVec[lumpVec[i].li_I];
      fVec[lumpVec[i].li_V2]  -= solVec[lumpVec[i].li_I];
      fVec[lumpVec[i].li_I ]  -= (solVec[lumpVec[i].li_V1]-solVec[lumpVec[i].li_V2]);

      // resistor
      double i0_res = (solVec[lumpVec[i].li_V2]-solVec[lumpVec[i].li_V3])*G;

      fVec[lumpVec[i].li_V2]  += i0_res;
      fVec[lumpVec[i].li_V3]  -= i0_res;
    }
  }
  else if ( model_.specialCase == TRANS_MOD_LC)
  {
    for (int i=0; i<numLumps; i++)
    {
      // inductor
      fVec[lumpVec[i].li_V1]  += solVec[lumpVec[i].li_I];
      fVec[lumpVec[i].li_V2]  -= solVec[lumpVec[i].li_I];
      fVec[lumpVec[i].li_I ]  -= (solVec[lumpVec[i].li_V1]-solVec[lumpVec[i].li_V2]);
    }
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEdQdx
//
// Purpose       : Loads the dQdx-matrix contributions for a single
//                 TransLine instance.
//
// Special Notes : The "Q" std::vector is part of a standard DAE formalism in
//                 which the system of equations is represented as:
//
//                 f(x) = dQ(x)/dt + F(x) - B(t) = 0
//
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 9/17/2013
//-----------------------------------------------------------------------------
bool Instance::loadDAEdQdx ()
{
  Linear::Matrix & dQdx = *(extData.dQdxMatrixPtr);
  if ( model_.specialCase == TRANS_MOD_RLC)
  {
    for (int i=0; i<numLumps; i++)
    {
      // inductor
      dQdx[lumpVec[i].li_I][lumpVec[i].offset_ii_ii] += L;

      // capacitor
      dQdx[lumpVec[i].li_V3][lumpVec[i].offset_v3_v3] += C;
    }
  }
  else if ( model_.specialCase == TRANS_MOD_LC)
  {
    for (int i=0; i<numLumps; i++)
    {
      // inductor
      dQdx[lumpVec[i].li_I][lumpVec[i].offset_ii_ii] += L;

      // capacitor
      dQdx[lumpVec[i].li_V2][lumpVec[i].offset_v2_v2] += C;
    }
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEdFdx ()
//
// Purpose       : Loads the F-vector contributions for a single
//                 TransLine instance.
//
// Special Notes :
//
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 9/17/2013
//-----------------------------------------------------------------------------
bool Instance::loadDAEdFdx ()
{
  Linear::Matrix & dFdx = *(extData.dFdxMatrixPtr);

  if ( model_.specialCase == TRANS_MOD_RLC)
  {
    for (int i=0; i<numLumps; i++)
    {
      // inductor
      dFdx[lumpVec[i].li_V1][lumpVec[i].offset_v1_ii] += 1.0;
      dFdx[lumpVec[i].li_V2][lumpVec[i].offset_v2_ii] -= 1.0;
      dFdx[lumpVec[i].li_I ][lumpVec[i].offset_ii_v1] -= 1.0;
      dFdx[lumpVec[i].li_I ][lumpVec[i].offset_ii_v2] += 1.0;

      // resistor
      dFdx[lumpVec[i].li_V2][lumpVec[i].offset_v2_v2] += G;
      dFdx[lumpVec[i].li_V2][lumpVec[i].offset_v2_v3] -= G;
      dFdx[lumpVec[i].li_V3][lumpVec[i].offset_v3_v2] -= G;
      dFdx[lumpVec[i].li_V3][lumpVec[i].offset_v3_v3] += G;
    }
  }
  else if ( model_.specialCase == TRANS_MOD_LC)
  {
    for (int i=0; i<numLumps; i++)
    {
      // inductor
      dFdx[lumpVec[i].li_V1][lumpVec[i].offset_v1_ii] += 1.0;
      dFdx[lumpVec[i].li_V2][lumpVec[i].offset_v2_ii] -= 1.0;
      dFdx[lumpVec[i].li_I ][lumpVec[i].offset_ii_v1] -= 1.0;
      dFdx[lumpVec[i].li_I ][lumpVec[i].offset_ii_v2] += 1.0;
    }
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::varTypes
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 9/17/2013
//-----------------------------------------------------------------------------
void Instance::varTypes( std::vector<char> & varTypeVec )
{
}


// Class Model

//-----------------------------------------------------------------------------
// Function      : Model::processParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 9/17/2013
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
// Creator       : Eric Keiter, SNL
// Creation Date : 9/17/2013
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
// Creator       : Eric Keiter, SNL
// Creation Date : 9/17/2013
//-----------------------------------------------------------------------------

Model::Model(
  const Configuration & configuration,
  const ModelBlock &    MB,
  const FactoryBlock &  factory_block)
  : DeviceModel(MB, configuration.getModelParameters(), factory_block),
    elevNumber(2),
    resist(0.0),
    induct(0.0),
    conduct(0.0),
    capac(0.0),

    elevNumberGiven(false),
    resistGiven(false),
    inductGiven(false),
    conductGiven(false),
    capacGiven(false),
    specialCase(TRANS_MOD_RLC)
{

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
// Function      : Model::Model
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 9/17/2013
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
// Creator       : Eric Keiter, SNL
// Creation Date : 9/17/2013
//-----------------------------------------------------------------------------

std::ostream &Model::printOutInstances(std::ostream &os) const
{
  std::vector<Instance*>::const_iterator iter;
  std::vector<Instance*>::const_iterator first = instanceContainer.begin();
  std::vector<Instance*>::const_iterator last  = instanceContainer.end();

  int i,isize;

  isize = instanceContainer.size();
  os << std::endl;
  os << "Number of TransLine instances: " << isize << std::endl;
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

//-----------------------------------------------------------------------------
// Function      : Xyce::Device::TRA::Master::loadDAEMatrices
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 11/26/08
//-----------------------------------------------------------------------------
///
/// Load DAE matrices for all resistor instances, regardless of model
///
/// @param dFdx matrix of derivatives of F vector with respect to solution
/// @param dQdx matrix of derivatives of Q vector with respect to solution
///
/// @return true on success
///
/// @note Because the resistor device re-implements the base-class
/// Master::loadDAEMatrices, the Instance::loadDAEdFdx method is
/// never called.  This method replaces those, and does the same work
/// but inside a loop over all resistor instances.
///
/// @see Xyce::Device::TRA::Instance::loadDAEdFdx
///
/// @author Eric Keiter, SNL
/// @date   11/26/08

bool Master::loadDAEMatrices(Linear::Matrix & dFdx, Linear::Matrix & dQdx, int loadType)
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
    Instance & tranlinei = *(*it);

    if ( tranlinei.model_.specialCase == TRANS_MOD_RLC)
    {
      for (int i=0; i<tranlinei.numLumps; i++)
      {
        // inductor
        dFdx[tranlinei.lumpVec[i].li_V1][tranlinei.lumpVec[i].offset_v1_ii] += 1.0;
        dFdx[tranlinei.lumpVec[i].li_V2][tranlinei.lumpVec[i].offset_v2_ii] -= 1.0;
        dFdx[tranlinei.lumpVec[i].li_I ][tranlinei.lumpVec[i].offset_ii_v1] -= 1.0;
        dFdx[tranlinei.lumpVec[i].li_I ][tranlinei.lumpVec[i].offset_ii_v2] += 1.0;
        dQdx[tranlinei.lumpVec[i].li_I][tranlinei.lumpVec[i].offset_ii_ii] += tranlinei.L;

        // resistor
        dFdx[tranlinei.lumpVec[i].li_V2][tranlinei.lumpVec[i].offset_v2_v2] += tranlinei.G;
        dFdx[tranlinei.lumpVec[i].li_V2][tranlinei.lumpVec[i].offset_v2_v3] -= tranlinei.G;
        dFdx[tranlinei.lumpVec[i].li_V3][tranlinei.lumpVec[i].offset_v3_v2] -= tranlinei.G;
        dFdx[tranlinei.lumpVec[i].li_V3][tranlinei.lumpVec[i].offset_v3_v3] += tranlinei.G;

        // capacitor
        dQdx[tranlinei.lumpVec[i].li_V3][tranlinei.lumpVec[i].offset_v3_v3] += tranlinei.C;
      }
    }
    else if ( tranlinei.model_.specialCase == TRANS_MOD_LC)
    {
      for (int i=0; i<tranlinei.numLumps; i++)
      {
        // inductor
        dFdx[tranlinei.lumpVec[i].li_V1][tranlinei.lumpVec[i].offset_v1_ii] += 1.0;
        dFdx[tranlinei.lumpVec[i].li_V2][tranlinei.lumpVec[i].offset_v2_ii] -= 1.0;
        dFdx[tranlinei.lumpVec[i].li_I ][tranlinei.lumpVec[i].offset_ii_v1] -= 1.0;
        dFdx[tranlinei.lumpVec[i].li_I ][tranlinei.lumpVec[i].offset_ii_v2] += 1.0;
        dQdx[tranlinei.lumpVec[i].li_I][tranlinei.lumpVec[i].offset_ii_ii] += tranlinei.L;

        // capacitor
        dQdx[tranlinei.lumpVec[i].li_V2][tranlinei.lumpVec[i].offset_v2_v2] += tranlinei.C;
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
// Creator       : Heidi Thornquist, SNL
// Creation Date : 11/26/08
//-----------------------------------------------------------------------------
bool Master::loadDAEVectors (double * solVec, double * fVec, double *qVec,  double * bVec, 
                             double * leadF, double * leadQ, double * junctionV, int loadType)
{
  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << " ----------------------------------" << std::endl;
    Xyce::dout() << " Master::loadDAEVectors: " << std::endl;
  }

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
    Instance & tranlinei = *(*it);

    // Load F vector
    if ( tranlinei.model_.specialCase == TRANS_MOD_RLC)
    {
      for (int i=0; i<tranlinei.numLumps; i++)
      {
        // inductor
        fVec[tranlinei.lumpVec[i].li_V1]  += solVec[tranlinei.lumpVec[i].li_I];
        fVec[tranlinei.lumpVec[i].li_V2]  -= solVec[tranlinei.lumpVec[i].li_I];
        fVec[tranlinei.lumpVec[i].li_I ]  -= (solVec[tranlinei.lumpVec[i].li_V1]-solVec[tranlinei.lumpVec[i].li_V2]);

        // resistor
        double i0_res = (solVec[tranlinei.lumpVec[i].li_V2]-solVec[tranlinei.lumpVec[i].li_V3])*tranlinei.G;
        fVec[tranlinei.lumpVec[i].li_V2]  += i0_res;
        fVec[tranlinei.lumpVec[i].li_V3]  -= i0_res;

        qVec[tranlinei.lumpVec[i].li_I] += tranlinei.L*solVec[tranlinei.lumpVec[i].li_I];
        qVec[tranlinei.lumpVec[i].li_V3] += tranlinei.C*solVec[tranlinei.lumpVec[i].li_V3];
      }
    }
    else if ( tranlinei.model_.specialCase == TRANS_MOD_LC)
    {
      for (int i=0; i<tranlinei.numLumps; i++)
      {
        // inductor
        fVec[tranlinei.lumpVec[i].li_V1]  += solVec[tranlinei.lumpVec[i].li_I];
        fVec[tranlinei.lumpVec[i].li_V2]  -= solVec[tranlinei.lumpVec[i].li_I];
        fVec[tranlinei.lumpVec[i].li_I ]  -= (solVec[tranlinei.lumpVec[i].li_V1]-solVec[tranlinei.lumpVec[i].li_V2]);

        qVec[tranlinei.lumpVec[i].li_I] += tranlinei.L*solVec[tranlinei.lumpVec[i].li_I];
        qVec[tranlinei.lumpVec[i].li_V2] += tranlinei.C*solVec[tranlinei.lumpVec[i].li_V2];
      }
    }
  }

  return true;
}

Device *Traits::factory(const Configuration &configuration, const FactoryBlock &factory_block)
{
  return new Master( configuration, factory_block, factory_block.solverState_, factory_block.deviceOptions_);
}

void
registerDevice(const DeviceCountMap& deviceMap, const std::set<int>& levelSet)
{
  if (deviceMap.empty() || (deviceMap.find("TRANSLINE")!=deviceMap.end()))
  {
    Config<Traits>::addConfiguration()
      .registerDevice("transline", 1)
      .registerModelType("transline", 1);
  }
}

} // namespace TransLine
} // namespace Device
} // namespace Xyce

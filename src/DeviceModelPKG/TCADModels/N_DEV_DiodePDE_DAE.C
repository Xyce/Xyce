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


//-------------------------------------------------------------------------
//
// Purpose        : One dimensional PDE device, new-DAE functions.
//
// Special Notes  :
//
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 05/13/05
//
//
//
//
//-------------------------------------------------------------------------

#include <Xyce_config.h>


// ----------  Standard Includes ----------
#include <iostream>
#include <N_UTL_Math.h>
#include <cstdio>

// ----------   Xyce Includes   ----------
#include <N_DEV_DiodePDE.h>
#include <N_DEV_ExternData.h>
#include <N_DEV_SolverState.h>
#include <N_DEV_DeviceOptions.h>
#include <N_DEV_MatrixLoadData.h>

#include <N_LAS_Vector.h>
#include <N_LAS_Matrix.h>
#include <N_LAS_System.h>
#include <N_LAS_Builder.h>

#include <N_NLS_TwoLevelPrintJac.h>

#include <N_ERH_ErrorMgr.h>
#include <N_DEV_Message.h>
// default number of mesh points:
#define NUM_MESH_POINTS 11

// default maximum number of nonzero entries in a matrix row
#define MAX_COLS_PER_ROW 10

namespace Xyce {
namespace Device {
namespace DiodePDE {

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEFVector
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/18/05
//-----------------------------------------------------------------------------
bool Instance::loadDAEFVector ()
{
  bool bsuccess = true;
  bool bs1;

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << "loadDAEFVector:  doubleDCOPStep="<<getSolverState().doubleDCOPStep;
    if (getSolverState().dcopFlag) Xyce::dout() << "  DCOP load" << std::endl;
    else                   Xyce::dout() << "  Transient load" << std::endl;
  }

  if ((getSolverState().dcopFlag) && getSolverState().doubleDCOPStep == 0)
  {
    equationSet = 0;
    bs1 = loadDAEFNonlinPoisson ();
  }
  else
  {
    equationSet = 1;

    if (getSolverState().twoLevelNewtonCouplingMode==Nonlinear::INNER_PROBLEM ||
        getSolverState().twoLevelNewtonCouplingMode==Nonlinear::FULL_PROBLEM)
    {
      bs1 = loadDAEFDDFormulation ();
    }
    else if (getSolverState().twoLevelNewtonCouplingMode==Nonlinear::OUTER_PROBLEM)
    {
      bs1 = loadDAEFExtractedConductance ();
    }
    else
    {
      DevelFatal(*this).in("Instance::loadDAEFVector")
        << "Invalid coupling Mode";
    }
  }

  return (bsuccess && bs1);
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEFNonlinPoisson
// Purpose       : Loads the F-vector the nonlinear poisson calculation.
//
// Special Notes : This should be identical to the original loadRHS function,
//                 for the nonlinear poisson.  All of that goes into the
//                 FVector, but with the opposite sign.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/18/05
//-----------------------------------------------------------------------------
bool Instance::loadDAEFNonlinPoisson ()
{
  double * fvec = extData.daeFVectorRawPtr;

  bool bsuccess = loadVecNLPoisson ( fvec );

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEFDDFormulation
// Purpose       : This function should be called from the loadDAEFVector
//                 function when solving the drift-diffusion equations.
// Special Notes :
//
// Scope         : private
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/18/05
//-----------------------------------------------------------------------------
bool Instance::loadDAEFDDFormulation ()
{
  double * fvec = extData.daeFVectorRawPtr;

  bool bsuccess = loadVecDDForm ( fvec );

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEFExtractedConductance
// Purpose       :
// Special Notes :
// Scope         : private
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 05/18/05
//-----------------------------------------------------------------------------
bool Instance::loadDAEFExtractedConductance ()
{
  Linear::Vector & fVec = *(extData.daeFVectorPtr) ;

  // KCL equations for the various connecting terminals:
  int iRow = 0;
  for (int iBC=0;iBC<bcVec.size();++iBC, ++iRow)
  {
    double voltLimFac = 0.0;
    if (getDeviceOptions().voltageLimiterFlag && voltLimFlag)
    {
      for (int iCol=0; iCol < numElectrodes ; ++iCol)
      {
        double vdiff = bcVec[iCol].Vckt_final - bcVec[iCol].Vckt_orig;
        vdiff *= scalingVars.V0;
        voltLimFac += vdiff * condVec[iRow][iCol];
      }
    }

    fVec[bcVec[iBC].lid] += bcVec[iBC].currentSum - voltLimFac;
  }

  for (int i=0;i<NX;++i)
  {
    fVec[ li_Vrowarray[i] ] = 0.0;
    fVec[ li_Nrowarray[i] ] = 0.0;
    fVec[ li_Prowarray[i] ] = 0.0;
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEQVector
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/18/05
//-----------------------------------------------------------------------------
bool Instance::loadDAEQVector ()
{
  bool bsuccess = true;
  bool bs1;

  Linear::Vector & qvec = *(extData.daeQVectorPtr);

  if ((getSolverState().dcopFlag) && getSolverState().doubleDCOPStep == 0)
  {
    equationSet = 0;

    bs1 = true; // no-op
  }
  else
  {
    equationSet = 1;

    if (getSolverState().twoLevelNewtonCouplingMode==Nonlinear::INNER_PROBLEM ||
        getSolverState().twoLevelNewtonCouplingMode==Nonlinear::FULL_PROBLEM)
    {
      bs1 = loadDAEQDDFormulation ();
    }
    else if (getSolverState().twoLevelNewtonCouplingMode==Nonlinear::OUTER_PROBLEM)
    {
      bs1 = loadDAEQExtractedConductance ();
    }
    else
    {
      DevelFatal(*this).in("Instance::loadDAEQVector")
        << "Invalid coupling Mode";
    }
  }

  return (bsuccess && bs1);
}


//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEQDDFormulation
// Purpose       :
// Special Notes :
// Scope         : private
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 05/18/05
//-----------------------------------------------------------------------------
bool Instance::loadDAEQDDFormulation ()
{
  Linear::Vector & qvec = *(extData.daeQVectorPtr);
  int i;

  // mesh points for the PDE problem:
  // only do the interior mesh points - boundaries have no dndt terms.
  for (i=1;i<NX-1;++i)
  {
    qvec[ li_Nrowarray[i] ] = -nnVec[i]*scalingVars.t0;
    qvec[ li_Prowarray[i] ] = -npVec[i]*scalingVars.t0;
  } // ip_iter, row loop...

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEQExtractedConductance
// Purpose       :
// Special Notes :
// Scope         : private
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 05/18/05
//-----------------------------------------------------------------------------
bool Instance::loadDAEQExtractedConductance ()
{
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEdFdx
// Purpose       : This function performs an analytic Jacobian matrix load for
//                 the diode-pde class, for the case of solving a nonlinear
//                 poisson equation.
// Special Notes :
//
// Scope         : private
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 05/18/05
//-----------------------------------------------------------------------------
bool Instance::loadDAEdFdx ()
{
  bool bsuccess;

  if ((getSolverState().dcopFlag) && getSolverState().doubleDCOPStep == 0)
  {
    bsuccess = loadDAEdFdxNonlinPoisson ();
  }
  else
  {
    if (getSolverState().twoLevelNewtonCouplingMode==Nonlinear::INNER_PROBLEM ||
        getSolverState().twoLevelNewtonCouplingMode==Nonlinear::FULL_PROBLEM)
    {
      bsuccess = loadDAEdFdxDDFormulation ();
    }
    else if (getSolverState().twoLevelNewtonCouplingMode==Nonlinear::OUTER_PROBLEM)
    {
      bsuccess = loadDAEdFdxExtractedConductance ();
    }
    else
    {
      Report::DevelFatal().in("Instance::loadDAEdFdx") << "Invalid coupling Mode" << numElectrodes;
    }
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEdFdxNonlinPoisson
// Purpose       : This function performs an analytic Jacobian matrix load for
//                 the diode-pde class, for the case of solving a nonlinear
//                 poisson equation.
// Special Notes :
//
// Scope         : private
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 05/18/05
//-----------------------------------------------------------------------------
bool Instance::loadDAEdFdxNonlinPoisson ()
{
  return loadMatNLPoisson( *(extData.dFdxMatrixPtr) );
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEdFdxDDFormulation
// Purpose       :
// Special Notes :
// Scope         : private
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 05/18/05
//-----------------------------------------------------------------------------
bool Instance::loadDAEdFdxDDFormulation ()
{
  return loadMatDDForm ( *(extData.dFdxMatrixPtr));
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEdFdxExtractedConductance
// Purpose       :
// Special Notes :
// Scope         : private
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 05/18/05
//-----------------------------------------------------------------------------
bool Instance::loadDAEdFdxExtractedConductance ()
{
  bool bsuccess = true;
  bool bs1 = true;

  Linear::Matrix & dFdx = *(extData.dFdxMatrixPtr);

  // first put 1's on the diagonals of all the mesh-rows:
  for (int i=0;i<NX;++i)
  {
    int Vrow = li_Vrowarray[i], Nrow = li_Nrowarray[i], Prow = li_Prowarray[i];
    std::vector<int> & Voff = li_Vcolarray[i];
    std::vector<int> & Noff = li_Ncolarray[i];
    std::vector<int> & Poff = li_Pcolarray[i];
    dFdx[Vrow][Voff[1]] =  1.0;
    dFdx[Nrow][Noff[1]] =  1.0;
    dFdx[Prow][Poff[1]] =  1.0;
  }

  // load the equivalent conductances matrix
  // the jacLID load is ordered confusingly, as the first entry in the jacStamp 
  // was the diagonal, and then all the off-diagonals follow.  
  for (int iRow=0; iRow < numElectrodes ; ++iRow)
  {
    int crossIndex=0;
    for (int iCol=0; iCol < numElectrodes ; ++iCol)
    {
      int iRowLID = bcVec[iRow].lid;
      int iColLID = bcVec[iRow].lidOffset;
      if (iCol!=iRow) { iColLID = bcVec[iRow].crossOffsets[crossIndex++]; }
      dFdx[iRowLID][iColLID] += condVec[iRow][iCol];
    }
  }

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
  {
    std::vector<std::string> names;
    for (int iE1 = 0; iE1 < numElectrodes; ++iE1) { names.push_back(bcVec[iE1].eName); }

    Xyce::Nonlinear::printJacobian(Xyce::dout(), outputName, names,condVec);
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEdQdx
// Purpose       :
// Special Notes :
// Scope         : private
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 05/18/05
//-----------------------------------------------------------------------------
bool Instance::loadDAEdQdx ()
{
  bool bsuccess;
  Linear::Matrix & dQdxMat = *(extData.dQdxMatrixPtr);

  if ((getSolverState().dcopFlag) && getSolverState().doubleDCOPStep == 0)
  {
    bsuccess = true; // no-op
  }
  else
  {
    if (getSolverState().twoLevelNewtonCouplingMode==Nonlinear::INNER_PROBLEM ||
        getSolverState().twoLevelNewtonCouplingMode==Nonlinear::FULL_PROBLEM)
    {
      bsuccess = loadDAEdQdxDDFormulation ();
    }
    else if (getSolverState().twoLevelNewtonCouplingMode==Nonlinear::OUTER_PROBLEM)
    {
      bsuccess = loadDAEdQdxExtractedConductance ();
    }
    else
    {
      Report::DevelFatal().in("Instance::loadDAEdQdx") << "Invalid coupling Mode " << numElectrodes;
    }
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEdQdxDDFormulation
// Purpose       :
// Special Notes :
// Scope         : private
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 05/18/05
//-----------------------------------------------------------------------------
bool Instance::loadDAEdQdxDDFormulation ()
{
  bool bsuccess = true;

  Linear::Matrix & dQdxMat = *(extData.dQdxMatrixPtr);

  // load the rows associated with the PDE mesh:
  // use direct matrix access:
  for (int i=1;i<NX-1;++i)
  {
    int li_Nrow = li_Nrowarray[i];
    int li_Prow = li_Prowarray[i];

    // Note::these terms need to be the SAME sign as the terms
    // in loadDAEQVector!

    // electron continuity row:
    // derivative w.r.t. nnVec[i  ]:
    dQdxMat[li_Nrow][li_Ncolarray[i][1]] = -scalingVars.t0;

    // hole continuity row:
    // derivative w.r.t. npVec[i  ]:
    dQdxMat[li_Prow][li_Pcolarray[i][1]] = -scalingVars.t0;
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEdQdxExtractedConductance
// Purpose       :
// Special Notes :
// Scope         : private
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 05/18/05
//-----------------------------------------------------------------------------
bool Instance::loadDAEdQdxExtractedConductance ()
{
  bool bsuccess = true;
  return bsuccess;
}

} // namespace DiodePDE
} // namespace Device
} // namespace Xyce

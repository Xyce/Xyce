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
// Purpose        : This file contains a lot of the
//                  implementation of the instance class for the two
//                  dimensional PDE based semiconductor device.
//
//                  Functions pertaining to the initial setup are in other
//                  files, as are functions relating to mesh handling and
//                  parameter handling.
//
// Special Notes  :
//
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 07/05/03
//
//
//
//
//-------------------------------------------------------------------------

#include <Xyce_config.h>

// ----------  Standard Includes ----------
#include <iostream>

// ----------   Xyce Includes   ----------
#include <N_DEV_2DPDE.h>
#include <N_DEV_SolverState.h>
#include <N_DEV_ExternData.h>
#include <N_DEV_DeviceOptions.h>
#include <N_DEV_MatrixLoadData.h>
#include <N_DEV_PDE_2DMesh.h>
#include <N_DEV_SourceData.h>

#include <N_LAS_Vector.h>
#include <N_LAS_Matrix.h>
#include <N_LAS_System.h>
#include <N_LAS_Builder.h>

#include <N_NLS_TwoLevelPrintJac.h>

#include <N_ERH_ErrorMgr.h>
#include <N_DEV_Message.h>
namespace Xyce {
namespace Device {
namespace TwoDPDE {

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEFVector
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/21/05
//-----------------------------------------------------------------------------
bool Instance::loadDAEFVector ()
{
  bool bsuccess = true;
  bool bs1;

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << "loadDAEFVector:  doubleDCOPStep="<<getSolverState().doubleDCOPStep;
    if (getSolverState().dcopFlag) 
    {
      Xyce::dout() << "  DCOP load" <<std::endl;
    }
    else                   
    {
      Xyce::dout() << "  Transient load" <<std::endl;
    }
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
// Creation Date : 6/21/05
//-----------------------------------------------------------------------------
bool Instance::loadDAEFNonlinPoisson ()
{
  return loadVecNLPoisson ( -1.0, extData.daeFVectorPtr );
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEFDDFormulation
// Purpose       : This function should be called from the loadDAEFVector
//                 function when solving the drift-diffusion equations.
// Special Notes :
//
// Scope         : private
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/21/05
//-----------------------------------------------------------------------------
bool Instance::loadDAEFDDFormulation ()
{
  // calcTerminalCurrents needs to be here b/c it is otherwise
  // called from updateSecondaryState
  calcTerminalCurrents ();
  return loadVecDDForm (-1.0,0.0, extData.daeFVectorPtr);
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEFExtractedConductance
// Purpose       :
// Special Notes :
// Scope         : private
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/21/05
//-----------------------------------------------------------------------------
bool Instance::loadDAEFExtractedConductance ()
{
  bool bsuccess = true;
  bool bs1 = true;
  int Vrow, Nrow, Prow;
  double coef;

  // KCL equations for the various connecting terminals:
  std::vector<DeviceInterfaceNode>::iterator firstDI = dIVec.begin();
  std::vector<DeviceInterfaceNode>::iterator lastDI  = dIVec.end ();
  std::vector<DeviceInterfaceNode>::iterator iterDI  = firstDI;

  int iRow = 0;
  for (;iterDI!=lastDI;++iterDI, ++iRow)
  {
    if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
    {
      Xyce::dout() << "name = " << iterDI->eName;
      Xyce::dout() << "  iRow = " << iRow << std::endl;
    }

    coef = iterDI->currentSum; 

    double voltLimFac = 0.0;
    if (getDeviceOptions().voltageLimiterFlag && voltLimFlag)
    {
      int iCol;
      for (iCol=0; iCol < numElectrodes ; ++iCol)
      {
        //if (dIVec[iCol].gid == -1) continue;

        double vdiff = dIVec[iCol].Vckt_final - dIVec[iCol].Vckt_orig;
        vdiff *= scalingVars.V0;
        voltLimFac += vdiff * condVec[iRow][iCol];

        if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
        {
          Xyce::dout() << "iCol = " << iCol;
          Xyce::dout() << "  vdiff = " << vdiff;
          Xyce::dout() << "  cond  = " << condVec[iRow][iCol];
          Xyce::dout() << "  voltLimFac = " << voltLimFac << std::endl;
        }
      }
    }

    (*extData.daeFVectorPtr)[iterDI->lid] += coef - voltLimFac;

    if (DEBUG_DEVICE)
    {
      (*extData.JdxpVectorPtr)[iterDI->lid] += voltLimFac;
    }

    if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
    {
      Xyce::dout() << "KCL for "<< iterDI->eName << ":\n";
      Xyce::dout() << " row       = " << iterDI->lid << "\n";
      Xyce::dout() << "currentSum = " << coef << "\n";
      Xyce::dout() << "voltLimFac = " << voltLimFac << std::endl;
    }
  } // end of DI loop

  // Load in the zeros...
  for (int i=0;i<numMeshPoints;++i)
  {
    if (boundarySten[i]) continue;

    Vrow = li_Vrowarray[i];
    Nrow = li_Nrowarray[i];
    Prow = li_Prowarray[i];
    (*extData.daeFVectorPtr)[Vrow] = 0.0;
    (*extData.daeFVectorPtr)[Nrow] = 0.0;
    (*extData.daeFVectorPtr)[Prow] = 0.0;
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEQVector
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/21/05
//-----------------------------------------------------------------------------
bool Instance::loadDAEQVector ()
{
  bool bsuccess = true;
  bool bs1;

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << "loadDAEQVector:  doubleDCOPStep="<<getSolverState().doubleDCOPStep;
    if (getSolverState().dcopFlag) Xyce::dout() << "  DCOP load" <<std::endl;
    else                   Xyce::dout() << "  Transient load" <<std::endl;
  }

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
// Creation Date : 6/21/05
//-----------------------------------------------------------------------------
bool Instance::loadDAEQDDFormulation ()
{
  bool bsuccess = true;
  bool bs1 = true;

  Linear::Vector * vecPtr = extData.daeQVectorPtr;

  int i;
  int Nrow, Prow;

  // mesh points for the PDE problem:
  for (i=0;i<numMeshPoints;++i)
  {
    if (boundarySten[i]) continue;

    Nrow = li_Nrowarray[i];
    Prow = li_Prowarray[i];
    (*vecPtr)[Nrow] = -nnVec[i]*scalingVars.t0;
    (*vecPtr)[Prow] = -npVec[i]*scalingVars.t0;
  }

  return bsuccess;
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
// Creation Date : 6/21/05
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
      DevelFatal(*this).in("Instance::loadDAEdFdx") << "Invalid coupling Mode" << numElectrodes;
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
// Creation Date : 6/21/05
//-----------------------------------------------------------------------------
bool Instance::loadDAEdFdxNonlinPoisson ()
{
  return loadMatNLPoisson( extData.dFdxMatrixPtr );
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEdFdxDDFormulation
// Purpose       :
// Special Notes :
// Scope         : private
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/21/05
//-----------------------------------------------------------------------------
bool Instance::loadDAEdFdxDDFormulation ()
{
  bool bsuccess = true;
  bool bs1 = true;
  // set up some of the partial derivative arrays:
  bs1 = pdRecombination ();   bsuccess = bsuccess && bs1;
  bs1 = pdElectronCurrent (); bsuccess = bsuccess && bs1;
  bs1 = pdHoleCurrent ();     bsuccess = bsuccess && bs1;
  bs1 = pdTerminalCurrents ();bsuccess = bsuccess && bs1;

  if ( !(getSolverState().twoLevelNewtonCouplingMode==Nonlinear::INNER_PROBLEM))
  {
    bs1 = loadMatKCLDDForm ( extData.dFdxMatrixPtr );
    bsuccess = bsuccess && bs1;
  }
  else
  {
    bs1 = loadMatCktTrivial ( extData.dFdxMatrixPtr );
    bsuccess = bsuccess && bs1;
  }

  bs1 = loadMatDDForm ( 0.0, extData.dFdxMatrixPtr );
  bsuccess = bsuccess && bs1;

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEdFdxExtractedConductance
// Purpose       :
// Special Notes :
// Scope         : private
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/21/05
//-----------------------------------------------------------------------------
bool Instance::loadDAEdFdxExtractedConductance ()
{
  Linear::Matrix & dFdx = *(extData.dFdxMatrixPtr);

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << std::endl;
    Xyce::dout() << section_divider << std::endl;
    Xyce::dout() << "Instance::loadDAEdFdxExtractedConductance" << std::endl;
  }

  // put 1's on the diagonals of all the mesh-rows:
  for (int i=0;i<numMeshPoints;++i)
  {
    if (boundarySten[i]) continue;

    int Vrow = li_Vrowarray[i], Nrow = li_Nrowarray[i], Prow = li_Prowarray[i];
    std::vector<int> & Voff = li_VoffsetArray[i];
    std::vector<int> & Noff = li_NoffsetArray[i];
    std::vector<int> & Poff = li_PoffsetArray[i];
    dFdx[Vrow][Voff[0]] =  1.0;
    dFdx[Nrow][Noff[0]] =  1.0;
    dFdx[Prow][Poff[0]] =  1.0;
  }

  // load the equivalent conductances matrix
  // the jacLID load is ordered confusingly, as the first entry in the jacStamp 
  // was the diagonal, and then all the off-diagonals follow.  
  for (int iRow=0; iRow < numElectrodes ; ++iRow)
  {
    int crossIndex=0;
    for (int iCol=0; iCol < numElectrodes ; ++iCol)
    {
      int iRowLID = dIVec[iRow].lid;
      int iColLID = dIVec[iRow].lidOffset;
      if (iCol!=iRow) { iColLID = dIVec[iRow].crossOffsets[crossIndex++]; }
      dFdx[iRowLID][iColLID] += condVec[iRow][iCol];
    }
  }

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
  {
    std::vector<std::string> names;
    for (int iE1 = 0; iE1 < numElectrodes; ++iE1) { names.push_back(dIVec[iE1].eName); }

    Xyce::Nonlinear::printJacobian(Xyce::dout(),outputName,names,condVec);
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEdQdx
// Purpose       :
// Special Notes :
// Scope         : private
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/21/05
//-----------------------------------------------------------------------------
bool Instance::loadDAEdQdx ()
{
  bool bsuccess;

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
      Report::DevelFatal().in("Instance::loadDAEdQdx") << "Invalid coupling Mode" << numElectrodes;
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
// Creation Date : 6/21/05
//-----------------------------------------------------------------------------
bool Instance::loadDAEdQdxDDFormulation ()
{
  bool bsuccess = true;
  bool bs1 = true;
  int i,j,count;
  int Nrow, Prow;

  int numCol = cols.size();
  if (vals.size () < cols.size()) numCol = vals.size();

  Linear::Matrix & QMatrix = (*extData.dQdxMatrixPtr);

  // load the rows associated with the PDE mesh:
  for (i=0;i<numMeshPoints;++i)
  {
    if (boundarySten[i]) continue;

    Nrow = li_Nrowarray[i];
    Prow = li_Prowarray[i];

    std::vector<int> & Noff = li_NoffsetArray[i];
    std::vector<int> & Poff = li_PoffsetArray[i];

    QMatrix[Nrow][Noff[0]] += -scalingVars.t0;
    QMatrix[Prow][Poff[0]] += -scalingVars.t0;

  } // mesh loop.

  return bsuccess;
}

} // namespace TwoDPDE
} // namespace Device
} // namespace Xyce


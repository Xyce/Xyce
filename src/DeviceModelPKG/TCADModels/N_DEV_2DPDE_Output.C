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

//-------------------------------------------------------------------------
//
// Purpose        : This file contains just the output functions of the
//                  N_DEV_2DPDE and N_DEV_2DPDEInstance classes.
//
// Special Notes  :
//
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 06/17/03
//
//
//
//
//-------------------------------------------------------------------------

#include <Xyce_config.h>

// ----------  Standard Includes ----------
#include <iostream>
#include <N_UTL_Math.h>

// ----------   Xyce Includes   ----------
#include <N_DEV_2DPDE.h>
#include <N_DEV_SolverState.h>
#include <N_DEV_DeviceOptions.h>
#include <N_DEV_Message.h>

#include <N_DEV_PDE_2DMesh.h>
#include <N_DEV_SGF_Interface.h>

#include <N_LAS_Vector.h>
#include <N_LAS_System.h>

namespace Xyce {
namespace Device {
namespace TwoDPDE {

//-----------------------------------------------------------------------------
// Function      : Instance::outputPlotFiles
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 07/22/03
//-----------------------------------------------------------------------------
bool Instance::outputPlotFiles(bool force_final_output)
{
  bool bsuccess = true;
  bool bs1 = true;
  bool skipOutput = false;

  // usually, don't bother outputting nonlinear Poisson result.
  if (equationSet == 0 && !outputNLPoisson)  return bsuccess;

  // If using output interval, check if enough time has passed to do
  // another output.  (only applies for transient - not DCOP).
  if ( !(getSolverState().dcopFlag) &&
  !(force_final_output) &&
    given("OUTPUTINTERVAL") )
  {
    double outMult = static_cast<double> (outputIndex);
    double nextOutputTime = outMult * outputInterval;

    if (nextOutputTime > getSolverState().currTime_)
    {
      skipOutput = true;
    }
  }

  // If this is a "forced" final output, make sure that it didn't already output.
  // This can happen if the output interval is an exact multiple of the
  // total simulation time.
  if (force_final_output && getSolverState().currTime_==lastOutputTime) skipOutput=true;

  if (skipOutput) return bsuccess;
  ++outputIndex;
  lastOutputTime = getSolverState().currTime_;

  if (DEBUG_DEVICE)
  {
    Xyce::dout() << std::endl << "Doing an output at time = " << getSolverState().currTime_ << std::endl;
  }

  if (tecplotLevel > 0)
  {
    bs1 = outputTecplot ();
    bsuccess = bsuccess && bs1;
  }

  if (tecplotLevel > 2)
  {
    bs1 = outputTecplotVectors ();
    bsuccess = bsuccess && bs1;
  }

  if (sgplotLevel > 0)
  {
    bs1 = outputSgplot ();
    bsuccess = bsuccess && bs1;
  }

  if (gnuplotLevel > 0)
  {
    bs1 = outputGnuplot ();
    bsuccess = bsuccess && bs1;
  }

  if (txtDataLevel > 0)
  {
    bs1 = outputTxtData ();
    bsuccess = bsuccess && bs1;
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::tecplotGeomOutput
//
// Purpose       : This function outputs the reactor geometry
//                 information for  a tecplot file.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/25/02
//-----------------------------------------------------------------------------
bool Instance::tecplotGeomOutput(FILE  *fp1)
{
  UINT iNodeA,iNodeB;

  // output the geometries:
  UINT iE;
  UINT iNumG = iNumPlotEdges/50 + 1;
  UINT iNumF = iNumPlotEdges - (iNumG-1)*50;

  UINT iGnum = 1;
  fprintf(fp1,"%s","\n  GEOMETRY M=GRID, C=BLACK, X= .00, Y= .00,");
  fprintf(fp1,"%s"," T=LINE, F=POINT, LT=0.8\n");

  if(iGnum == iNumG) fprintf(fp1,"\t%d\n",iNumF);
  else               fprintf(fp1,"\t%d\n",50);

  double x1,y1,x2,y2;
  UINT itmp = 2;
  UINT icount = 0;

  for(iE = 0;iE<numMeshEdges;++iE)
  {
    if(aiEdge[iE] == 1)
    {
      mEdge * edgePtr = meshContainerPtr->getEdge(iE);

      ++icount;
      iNodeA = edgePtr->inodeA;
      iNodeB = edgePtr->inodeB;
      x1 = xVec[iNodeA];
      y1 = yVec[iNodeA];
      x2 = xVec[iNodeB];
      y2 = yVec[iNodeB];

      if (variablesScaled)
      {
        x1 = x1 * scalingVars.x0;
        x2 = x2 * scalingVars.x0;
        y1 = y1 * scalingVars.x0;
        y2 = y2 * scalingVars.x0;
      }

      fprintf(fp1,"%4d\n%11.3e %11.3e\n%11.3e %11.3e\n",
		      itmp,x1,y1,x2,y2);
    }

    if(icount >= 50)
    {
      icount = 0;
      ++iGnum;
      if(iGnum == iNumG)
      {
        if(iNumF != 0)
        {
          fprintf(fp1,"%s","\n  GEOMETRY M=GRID, C=BLACK, X=    .00,");
          fprintf(fp1,"%s"," Y=    .00,");
          fprintf(fp1,"%s"," T=LINE, F=POINT, LT=0.8\n");
          fprintf(fp1,"\t%d\n",iNumF);}
        }
        else
        {
          fprintf(fp1,"%s","\n  GEOMETRY M=GRID, C=BLACK, X=    .00,");
          fprintf(fp1,"%s"," Y=    .00,");
          fprintf(fp1,"%s"," T=LINE, F=POINT, LT=0.8\n");
          fprintf(fp1,"\t%d\n",50);
        }
     }
  }

  fprintf(fp1,"%s","\n");

  // Now do the "noflux" edges.
  iNumG = iNumPlotEdges_nf/50 + 1;
  iNumF = iNumPlotEdges_nf - (iNumG-1)*50;

  iGnum = 1;
  fprintf(fp1,"%s","\n  GEOMETRY M=GRID, C=RED, X=    .00, Y=    .00,");
  fprintf(fp1,"%s"," T=LINE, F=POINT, LT=0.2\n");

  if(iGnum == iNumG) fprintf(fp1,"\t%d\n",iNumF);
  else               fprintf(fp1,"\t%d\n",50);

  itmp = 2;
  icount = 0;

  for(iE = 0;iE<numMeshEdges;++iE)
  {
    if(aiEdge_nf[iE] == 1)
    {
      mEdge * edgePtr = meshContainerPtr->getEdge(iE);

      ++icount;
      iNodeA = edgePtr->inodeA;
      iNodeB = edgePtr->inodeB;
      x1 = xVec[iNodeA];
      y1 = yVec[iNodeA];
      x2 = xVec[iNodeB];
      y2 = yVec[iNodeB];

      if (variablesScaled)
      {
        x1 = x1 * scalingVars.x0;
        x2 = x2 * scalingVars.x0;
        y1 = y1 * scalingVars.x0;
        y2 = y2 * scalingVars.x0;
      }

      fprintf(fp1,"%4d\n%11.3e %11.3e\n%11.3e %11.3e\n",
		      itmp,x1,y1,x2,y2);
    }

    if(icount >= 50)
    {
      icount = 0;
      ++iGnum;
      if(iGnum == iNumG)
      {
        if(iNumF != 0)
        {
          fprintf(fp1,"%s","\n  GEOMETRY M=GRID, C=RED, X=    .00,");
          fprintf(fp1,"%s"," Y=    .00,");
          fprintf(fp1,"%s"," T=LINE, F=POINT, LT=0.2\n");
          fprintf(fp1,"\t%d\n",iNumF);}
        }
        else
        {
          fprintf(fp1,"%s","\n  GEOMETRY M=GRID, C=RED, X=    .00,");
          fprintf(fp1,"%s"," Y=    .00,");
          fprintf(fp1,"%s"," T=LINE, F=POINT, LT=0.2\n");
          fprintf(fp1,"\t%d\n",50);
        }
     }
  }

  fprintf(fp1,"%s","\n");

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::outputTecplot
// Purpose       :
// Special Notes : If tecplot level is set to 1, then output each dataset
//                 in a separate file.  If not, then append to a single file.
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 11/14/01
//-----------------------------------------------------------------------------
bool Instance::outputTecplot ()
{
  bool bsuccess = true;
  bool bs1 = true;

  int i;
  char filename[32];   for(i=0;i<32;++i) filename[i] = static_cast<char>(0);

  if (tecplotLevel == 1)
  {
    sprintf(filename,"%s_%03d.dat",outputName.c_str(),callsOTEC);
  }
  else
  {
    sprintf(filename,"%s.dat",outputName.c_str());
  }

  double time = getSolverState().currTime_;

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << std::endl;
    Xyce::dout() << section_divider << std::endl;
    Xyce::dout() << "In Instance::outputTecplot.  filename = ";
    Xyce::dout() << std::string(filename);
    Xyce::dout() << std::endl;
  }

  FILE *fp1;

  if (tecplotLevel == 1)
  {
    fp1 = fopen(filename,"w");
  }
  else
  {
    if (callsOTEC <= 0)
    {
      fp1 = fopen(filename,"w");
    }
    else
    {
      fp1 = fopen(filename,"a");
    }
  }

  double Ex, Ey;
  double V2x, V2y;
  double V1x, V1y;
  double V0x, V0y;

  double dDx_dt, dDy_dt;
  double D2x, D2y;
  double D1x, D1y;
  double D0x, D0y;

  double xLoc, yLoc;
  double dx1, dy1;

  if (tecplotLevel == 1)
  {
    if (equationSet == 0)
    {
      fprintf(fp1,
    " TITLE = \"Spatially Dependent data for 2D PDE: %s  time = %12.4e seconds. equation set = nonlinear Poisson\",\n",
	outputName.c_str(),time);
    }
    else
    {
      fprintf(fp1,
    " TITLE = \"Spatially Dependent data for 2D PDE: %s  time = %12.4e seconds. equation set = drift diffusion\",\n",
	outputName.c_str(),time);
    }
  }
  else
  {
    if (callsOTEC <= 0)
    {
      fprintf(fp1,
       " TITLE = \"Spatially Dependent data for 2D PDE: %s  time = %12.4e seconds.\",\n",
      outputName.c_str(),time);
    }
  }

  if (callsOTEC <= 0 || tecplotLevel == 1)
  {
    if (cylGeomFlag)
      fprintf(fp1,"%s","\tVARIABLES = \"R (cm) \",\"Z (cm)\",\n");
    else
      fprintf(fp1,"%s","\tVARIABLES = \"X (cm) \",\"Y (cm)\",\n");

    fprintf(fp1,"%s","\t    \"V \",\n");
    fprintf(fp1,"%s","\t    \"nn (electron dens.) \",\n");
    fprintf(fp1,"%s","\t    \"np (hole dens.) \",\n");
    fprintf(fp1,"%s","\t    \"Dopant dens. \",\n");
    fprintf(fp1,"%s","\t    \"abs(Dopant dens.)\",\n");
    fprintf(fp1,"%s","\t    \"total density\",\n");
    fprintf(fp1,"%s","\t    \"electron lifetime \",\n");
    fprintf(fp1,"%s","\t    \"hole lifetime \",\n");
    fprintf(fp1,"%s","\t    \"electron mobility \",\n");
    fprintf(fp1,"%s","\t    \"hole mobility \",\n");

    if (tecplotLevel > 1)
    {
      fprintf(fp1,"%s","\t    \"Ex \",\n");
      fprintf(fp1,"%s","\t    \"Ey \",\n");
      fprintf(fp1,"%s","\t    \"Emag \",\n");
      fprintf(fp1,"%s","\t    \"dDx_dt \",\n");
      fprintf(fp1,"%s","\t    \"dDy_dt \",\n");
      fprintf(fp1,"%s","\t    \"dDmag_dt \",\n");
    }

    fprintf(fp1,"%s","\t    \"Recombination \",\n");
    fprintf(fp1,"%s","\t    \"photogen  \",\n");
    fprintf(fp1,"%s","\t    \"total src \",\n");
  }

  fprintf(fp1,"\tZONE F=FEPOINT,ET=QUADRILATERAL,N=%d,E=%d",
    numMeshPoints, numMeshCells);

  if (getSolverState().dcopFlag)
  {
    fprintf(fp1,"  T = \"DCOP step = %d\" \n", callsOTEC);
  }
  else
  {
    fprintf(fp1,"  T = \"time step = %d time = %12.4e seconds\" \n", callsOTEC , time);
  }

  if (variablesScaled)
  {
    for (i=0;i<numMeshPoints;++i)
    {
      xLoc = xVec[i];
      yLoc = yVec[i];

      double xLocS = xLoc * scalingVars.x0;
      double yLocS = yLoc * scalingVars.x0;

      fprintf(fp1,"  %12.4e",xLocS);
      fprintf(fp1,"  %12.4e",yLocS);
      fprintf(fp1,"  %12.4e",VVec[i]*scalingVars.V0);
      fprintf(fp1,"  %12.4e",nnVec[i]*scalingVars.C0);
      fprintf(fp1,"  %12.4e",npVec[i]*scalingVars.C0);
      double C = CVec[i]*scalingVars.C0;
      fprintf(fp1,"  %12.4e",C);
      fprintf(fp1,"  %12.4e",fabs(C));
      double totCharge = (npVec[i]*scalingVars.C0-nnVec[i]*scalingVars.C0+C);
      fprintf(fp1,"  %12.4e",totCharge);
      fprintf(fp1,"  %12.4e",tnVec[i]*scalingVars.t0);
      fprintf(fp1,"  %12.4e",tpVec[i]*scalingVars.t0);
      fprintf(fp1,"  %12.4e",unVec[i]*scalingVars.u0);
      fprintf(fp1,"  %12.4e",upVec[i]*scalingVars.u0);

      if (tecplotLevel > 1)
      {
        dx1 = dy1 = 0.5 * minDXVec[i];

        mLabel * lPtr = meshContainerPtr->getLabel(labelIndex[i]);
        if (lPtr->uType == TYPE_REGION)
        {
          // electric field:
          V2x = scalingVars.V0* (meshContainerPtr->interp(&(VVec[0]),xLocS+dx1, yLocS));
          V1x = scalingVars.V0*VVec[i];
          V0x = scalingVars.V0* (meshContainerPtr->interp(&(VVec[0]),xLocS-dx1, yLocS));
          Ex = -0.5 * ( (V2x-V1x)+(V1x-V0x) )/dx1;

          V2y = scalingVars.V0* (meshContainerPtr->interp(&(VVec[0]),xLocS, yLocS+dy1));
          V1y = scalingVars.V0*VVec[i];
          V0y = scalingVars.V0* (meshContainerPtr->interp(&(VVec[0]),xLocS, yLocS-dy1));
          Ey = -0.5 * ( (V2y-V1y)+(V1y-V0y) )/dy1;


          // displacement current:
          D2x = (meshContainerPtr->interp(&(displPotential[0]),xLocS+dx1, yLocS));
          D1x = displPotential[i];
          D0x = (meshContainerPtr->interp(&(displPotential[0]),xLocS-dx1, yLocS));

          dDx_dt = -0.5 * ( (D2x-D1x)+(D1x-D0x) )/dx1;
          dDx_dt *= eSi * e0;

          D2y = (meshContainerPtr->interp(&(displPotential[0]),xLocS, yLocS+dy1));
          D1y = displPotential[i];
          D0y = (meshContainerPtr->interp(&(displPotential[0]),xLocS, yLocS-dy1));

          dDy_dt = -0.5 * ( (D2y-D1y)+(D1y-D0y) )/dy1;
          dDy_dt *= eSi * e0;

        }
        else
        {
          Ex     = 0.0;
          Ey     = 0.0;
          dDx_dt = 0.0;
          dDy_dt = 0.0;
        }

        fprintf(fp1,"  %12.4e",Ex);
        fprintf(fp1,"  %12.4e",Ey);
        fprintf(fp1,"  %12.4e",sqrt(Ex*Ex + Ey*Ey) );

        fprintf(fp1,"  %12.4e",dDx_dt);
        fprintf(fp1,"  %12.4e",dDy_dt);
        fprintf(fp1,"  %12.4e",sqrt(dDx_dt*dDx_dt + dDy_dt*dDy_dt) );
      }

      fprintf(fp1,"  %12.4e",RVec[i]*scalingVars.R0);
      fprintf(fp1,"  %12.4e",SVec[i]*scalingVars.R0);
      fprintf(fp1,"  %12.4e",totSrcVec[i]*scalingVars.R0);
      fprintf(fp1,"%s","\n");
    }
  }
  else
  {
    for (i=0;i<numMeshPoints;++i)
    {
      xLoc = xVec[i];
      yLoc = yVec[i];

      fprintf(fp1,"  %12.4e",xLoc);
      fprintf(fp1,"  %12.4e",yLoc);
      fprintf(fp1,"  %12.4e",VVec[i]);
      fprintf(fp1,"  %12.4e",nnVec[i]);
      fprintf(fp1,"  %12.4e",npVec[i]);
      double C = CVec[i];
      fprintf(fp1,"  %12.4e",C);
      fprintf(fp1,"  %12.4e",fabs(C));
      double totCharge = (npVec[i]-nnVec[i]+C);
      fprintf(fp1,"  %12.4e",totCharge);
      fprintf(fp1,"  %12.4e",tnVec[i]);
      fprintf(fp1,"  %12.4e",tpVec[i]);
      fprintf(fp1,"  %12.4e",unVec[i]);
      fprintf(fp1,"  %12.4e",upVec[i]);

      if (tecplotLevel > 1)
      {
        dx1 = dy1 = 0.5 * minDXVec[i];
        mLabel * lPtr = meshContainerPtr->getLabel(labelIndex[i]);

        if (lPtr->uType == TYPE_REGION)
        {
          // electric field:
          V2x =  (meshContainerPtr->interp(&(VVec[0]),xLoc+dx1, yLoc));
          V1x = VVec[i];
          V0x =  (meshContainerPtr->interp(&(VVec[0]),xLoc-dx1, yLoc));
          Ex = -0.5 * ( (V2x-V1x)+(V1x-V0x) )/dx1;

          V2y =  (meshContainerPtr->interp(&(VVec[0]),xLoc, yLoc+dy1));
          V1y = VVec[i];
          V0y =  (meshContainerPtr->interp(&(VVec[0]),xLoc, yLoc-dy1));
          Ey = -0.5 * ( (V2y-V1y)+(V1y-V0y) )/dy1;

          // displacement current:
          D2x = (meshContainerPtr->interp(&(displPotential[0]),xLoc+dx1, yLoc));
          D1x = displPotential[i];
          D0x = (meshContainerPtr->interp(&(displPotential[0]),xLoc-dx1, yLoc));

          dDx_dt = -0.5 * ( (D2x-D1x)+(D1x-D0x) )/dx1;
          dDx_dt *= eSi * e0;

          D2y = (meshContainerPtr->interp(&(displPotential[0]),xLoc, yLoc+dy1));
          D1y = displPotential[i];
          D0y = (meshContainerPtr->interp(&(displPotential[0]),xLoc, yLoc-dy1));

          dDy_dt = -0.5 * ( (D2y-D1y)+(D1y-D0y) )/dy1;
          dDy_dt *= eSi * e0;
        }
        else
        {
          Ex     = 0.0;
          Ey     = 0.0;
          dDx_dt = 0.0;
          dDy_dt = 0.0;
        }

        fprintf(fp1,"  %12.4e",Ex);
        fprintf(fp1,"  %12.4e",Ey);
        fprintf(fp1,"  %12.4e",sqrt(Ex*Ex + Ey*Ey) );

        fprintf(fp1,"  %12.4e",dDx_dt);
        fprintf(fp1,"  %12.4e",dDy_dt);
        fprintf(fp1,"  %12.4e",sqrt(dDx_dt*dDx_dt + dDy_dt*dDy_dt) );
      }

      fprintf(fp1,"  %12.4e",RVec[i]);
      fprintf(fp1,"  %12.4e",SVec[i]);
      fprintf(fp1,"  %12.4e",totSrcVec[i]);
      fprintf(fp1,"%s","\n");
    }
  }

  fprintf(fp1,"%s","\n");

  for(UINT iTri=0; iTri<numMeshCells; ++iTri)
  {
    UINT inodeA,inodeB,inodeC,inodeD;

    mCell * cellPtr = meshContainerPtr->getCell(iTri);

    inodeA = cellPtr->inodeA;
    inodeB = cellPtr->inodeB;
    inodeC = cellPtr->inodeC;
    inodeD = cellPtr->inodeD;

    if(inodeD == -1u) inodeD = inodeC;

    if((inodeA == -1u || inodeB == -1u ||
        inodeC == -1u || inodeD == -1u) ||
       (inodeA == inodeB) || (inodeA == inodeC) ||
       (inodeB == inodeC))
    {
      fprintf(stdout,"%s","Error in ::outputTecplot\n");
      fprintf(stdout,"inodeA = %d\n",inodeA);
      fprintf(stdout,"inodeB = %d\n",inodeB);
      fprintf(stdout,"inodeC = %d\n",inodeC);
      fprintf(stdout,"inodeD = %d\n",inodeD);
    }

    // output
    fprintf(fp1,"%d %d %d %d\n",
      inodeA+1,
      inodeB+1,
      inodeC+1,
      inodeD+1);
  }

  if (tecplotLevel >= 1)
  {
    if (callsOTEC<=0)
    {
      bs1 = tecplotGeomOutput(fp1);
      bsuccess = bsuccess && bs1;
    }
  }

  ++callsOTEC;
  fclose(fp1);

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::outputTecplotVectors
// Purpose       : This function outputs a tecplot file which contains vector
//                 data.  It does this by performing interpolations on the
//                 unstructured grid.  As such, for testing purposed it also
//                 outputs some interpolated scalar data.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 11/14/01
//-----------------------------------------------------------------------------
bool Instance::outputTecplotVectors ()
{
  bool bsuccess = true;

  int i, j;
  char filename[32];   for(i=0;i<32;++i) filename[i] = static_cast<char>(0);

  sprintf(filename,"%s_%03dVec.dat",outputName.c_str(),callsOTECvec);
  ++callsOTECvec;

  double time = getSolverState().currTime_;

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << std::endl;
    Xyce::dout() << section_divider << std::endl;
    Xyce::dout() << "In Instance::outputTecplotVectors.  filename = ";
    Xyce::dout() << std::string(filename);
    Xyce::dout() << std::endl;
  }

  FILE *fp1 = fopen(filename,"w");

  // set up the cartesian grid with which to perform interpolations:
  double tmp;
  double xMax = meshContainerPtr->getXMax ();
  double xMin = meshContainerPtr->getXMin ();
  if (xMin > xMax)
  {
    tmp = xMax;
    xMax = xMin;
    xMin = tmp;
  }

  double yMax = meshContainerPtr->getYMax ();
  double yMin = meshContainerPtr->getYMin ();
  if (yMin > yMax)
  {
    tmp = yMax;
    yMax = yMin;
    yMin = tmp;
  }

  double DeltaX = xMax - xMin;
  double DeltaY = yMax - yMin;

  int inx;
  int iny;

  int iNumMeshPoints = interpGridSize;

  // Have one (the shorter one) of the axes have iNumMeshPoints "mesh" points.
  if (DeltaX < DeltaY)
  { inx = iNumMeshPoints; iny = static_cast<int> (iNumMeshPoints * fabs(DeltaY/DeltaX)); }
  else
  { iny = iNumMeshPoints; inx = static_cast<int> (iNumMeshPoints * fabs(DeltaX/DeltaY)); }

  double dx  = fabs(DeltaX/static_cast<double>(inx));
  double dx1 = dx*0.5;
  double dy  = fabs(DeltaY/static_cast<double>(iny));
  double dy1 = dy*0.5;

  double Ex, Ey;
  double V2x, V2y;
  double V1x, V1y;
  double V0x, V0y;

  double dDx_dt, dDy_dt;
  double D2x, D2y;
  double D1x, D1y;
  double D0x, D0y;

  if (DEBUG_DEVICE)
  {
    Xyce::dout() << "outputTecplotVectors:\n";
  Xyce::dout() << "DeltaX = " << DeltaX << std::endl;
  Xyce::dout() << "DeltaY = " << DeltaY << std::endl;

  Xyce::dout() << "xMax   = " << xMax   << std::endl;
  Xyce::dout() << "xMin   = " << xMin   << std::endl;
  Xyce::dout() << "yMax   = " << yMax   << std::endl;
  Xyce::dout() << "yMin   = " << yMin   << std::endl;

  Xyce::dout() << "dx     = " << dx  << std::endl;
  Xyce::dout() << "dy     = " << dy  << std::endl;

  Xyce::dout() << "inx    = " << inx << std::endl;
  Xyce::dout() << "iny    = " << iny << std::endl;
  }

  if (equationSet == 0)
  {
    fprintf(fp1,
  " TITLE = \"Spatially Dependent vector data for 2D PDE: %s  time = %12.4e seconds. equation set = nonlinear Poisson\",\n",
      outputName.c_str(),time);
  }
  else
  {
    fprintf(fp1,
  " TITLE = \"Spatially Dependent vector data for 2D PDE: %s  time = %12.4e seconds. equation set = drift diffusion\",\n",
      outputName.c_str(),time);
  }

  if (cylGeomFlag)
    fprintf(fp1,"%s","\tVARIABLES = \"R (cm) \",\"Z (cm)\",\n");
  else
    fprintf(fp1,"%s","\tVARIABLES = \"X (cm) \",\"Y (cm)\",\n");

  fprintf(fp1,"%s","\t    \"V \", \n");
  fprintf(fp1,"%s","\t    \"Ex \", \n");
  fprintf(fp1,"%s","\t    \"Ey \", \n");
  fprintf(fp1,"%s","\t    \"Emag \",\n");
  fprintf(fp1,"%s","\t    \"dDx_dt \", \n");
  fprintf(fp1,"%s","\t    \"dDy_dt \", \n");
  fprintf(fp1,"%s","\t    \"Dmag \",\n");

  fprintf(fp1,"\tZONE I=%d, J=%d, F=POINT\n", iny+1, inx+1);

  double xLoc, yLoc;

  if (variablesScaled)
  {
    for (i=0;i<inx+1;++i)
    {
      xLoc = xMin + static_cast<double>(i) * dx;

      for (j=0;j<iny+1;++j)
      {
        yLoc = yMin + static_cast<double>(j) * dy;

        fprintf(fp1,"  %12.4e",xLoc);
        fprintf(fp1,"  %12.4e",yLoc);
        double Vtmp = meshContainerPtr->interp(&(VVec[0]),xLoc,yLoc);
        double Dtmp = meshContainerPtr->interp(&(displPotential[0]),xLoc,yLoc);
        fprintf(fp1,"  %12.4e",Vtmp*scalingVars.V0);

        if (i>0 && i<inx)
        {
          // electric field
          V2x = scalingVars.V0* (meshContainerPtr->interp(&(VVec[0]),xLoc+dx1, yLoc));
          V1x = scalingVars.V0*Vtmp;
          V0x = scalingVars.V0* (meshContainerPtr->interp(&(VVec[0]),xLoc-dx1, yLoc));
          Ex = -0.5 * ( (V2x-V1x)+(V1x-V0x) )/dx1;

          // displacement current:
          D2x = (meshContainerPtr->interp(&(displPotential[0]),xLoc+dx1, yLoc));
          D1x = Dtmp;
          D0x = (meshContainerPtr->interp(&(displPotential[0]),xLoc-dx1, yLoc));

          dDx_dt = -0.5 * ( (D2x-D1x)+(D1x-D0x) )/dx1;
          dDx_dt *= eSi * e0;
        }
        else
        {
          Ex = 0.0;
          dDx_dt = 0.0;
        }

        if (j>0 && j<iny)
        {
          // electric field:
          V2y = scalingVars.V0* (meshContainerPtr->interp(&(VVec[0]),xLoc, yLoc+dy1));
          V1y = scalingVars.V0*Vtmp;
          V0y = scalingVars.V0* (meshContainerPtr->interp(&(VVec[0]),xLoc, yLoc-dy1));
          Ey = -0.5 * ( (V2y-V1y)+(V1y-V0y) )/dy1;

          // displacement current:
          D2y = (meshContainerPtr->interp(&(displPotential[0]),xLoc, yLoc+dy1));
          D1y = Dtmp;
          D0y = (meshContainerPtr->interp(&(displPotential[0]),xLoc, yLoc-dy1));

          dDy_dt = -0.5 * ( (D2y-D1y)+(D1y-D0y) )/dy1;
          dDy_dt *= eSi * e0;
        }
        else
        {
          Ey = 0.0;
          dDy_dt = 0.0;
        }

        fprintf(fp1,"  %12.4e",Ex);
        fprintf(fp1,"  %12.4e",Ey);
        fprintf(fp1,"  %12.4e",sqrt(Ex*Ex + Ey*Ey) );

        fprintf(fp1,"  %12.4e",dDx_dt);
        fprintf(fp1,"  %12.4e",dDy_dt);
        fprintf(fp1,"  %12.4e",sqrt(dDx_dt*dDx_dt + dDy_dt*dDy_dt) );
        fprintf(fp1,"%s","\n");
      }
    }
  }
  else
  {
    for (i=0;i<inx+1;++i)
    {
      xLoc = xMin + static_cast<double>(i) * dx;
      for (j=0;j<iny+1;++j)
      {
        yLoc = yMin + static_cast<double>(j) * dy;

        fprintf(fp1,"  %12.4e",xLoc);
        fprintf(fp1,"  %12.4e",yLoc);
        double Vtmp = meshContainerPtr->interp(&(VVec[0]),xLoc,yLoc);
        double Dtmp = meshContainerPtr->interp(&(displPotential[0]),xLoc,yLoc);
        fprintf(fp1,"  %12.4e",Vtmp);

        if (i>0 && i<inx)
        {
          // electric field:
          V2x = (meshContainerPtr->interp(&(VVec[0]),xLoc+dx1, yLoc));
          V1x = Vtmp;
          V0x = (meshContainerPtr->interp(&(VVec[0]),xLoc-dx1, yLoc));
          Ex = -0.5 * ( (V2x-V1x)+(V1x-V0x) )/dx1;

          // displacement current:
          D2x = (meshContainerPtr->interp(&(displPotential[0]),xLoc+dx1, yLoc));
          D1x = Dtmp;
          D0x = (meshContainerPtr->interp(&(displPotential[0]),xLoc-dx1, yLoc));

          dDx_dt = -0.5 * ( (D2x-D1x)+(D1x-D0x) )/dx1;
          dDx_dt *= eSi * e0;
        }
        else
        {
          Ex     = 0.0;
          dDx_dt = 0.0;
        }

        if (j>0 && j<iny)
        {
          // electric field:
          V2y = (meshContainerPtr->interp(&(VVec[0]),xLoc, yLoc+dy1));
          V1y = Vtmp;
          V0y = (meshContainerPtr->interp(&(VVec[0]),xLoc, yLoc-dy1));
          Ey = -0.5 * ( (V2y-V1y)+(V1y-V0y) )/dy1;

          // displacement current:
          D2y = (meshContainerPtr->interp(&(displPotential[0]),xLoc, yLoc+dy1));
          D1y = Dtmp;
          D0y = (meshContainerPtr->interp(&(displPotential[0]),xLoc, yLoc-dy1));

          dDy_dt = -0.5 * ( (D2y-D1y)+(D1y-D0y) )/dy1;
          dDy_dt *= eSi * e0;
        }
        else
        {
          Ey     = 0.0;
          dDy_dt = 0.0;
        }

        fprintf(fp1,"  %12.4e",Ex);
        fprintf(fp1,"  %12.4e",Ey);
        fprintf(fp1,"  %12.4e",sqrt(Ex*Ex + Ey*Ey) );

        fprintf(fp1,"  %12.4e",dDx_dt);
        fprintf(fp1,"  %12.4e",dDy_dt);
        fprintf(fp1,"  %12.4e",sqrt(dDx_dt*dDx_dt + dDy_dt*dDy_dt) );
        fprintf(fp1,"%s","\n");
      }
    }
  }

  // The geometry is based on the original mesh, which is output from "outputTecplot" as
  // 	ZONE F=FEPOINT,ET=QUADRILATERAL,N=#,E=#
  //
  // This vector output file is based on an interpolated cartesian mesh
  // as ZONE I=%d, J=%d, F=POINT
  //
  // So, the "tecplotGeomOutput" function as written doesn't work for the
  // vector data output file.
  //bool bs1 = true;
  //bs1 = tecplotGeomOutput(fp1);
  //bsuccess = bsuccess && bs1;

  fclose(fp1);

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::outputGnuplot
//
// Purpose       : This function outputs very simple column data that can
//                 be read by gnuplot, in the form x,y,V,n,p,C.
//
// Special Notes : It is easy to do this if we are running with
//                 the internal mesh, which is cartesian in appearance.
//
//                 If running without a cartesian mesh, then this function 
//                 isn't useful, as it would be neccessary to do interpolations.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/15/03
//-----------------------------------------------------------------------------
bool Instance::outputGnuplot ()
{
  bool bsuccess = true;

  int i, j;
  char filename[32];   for(i=0;i<32;++i) filename[i] = static_cast<char>(0);

  sprintf(filename,"%s_%03dGnu.dat",outputName.c_str(),callsOGNU);
  ++callsOGNU;

  if (!given("NX") || !given("NY"))
  {
    bsuccess = false;
    UserWarning(*this) << "Gnuplot only works if using the internal mesh.";
  }
  else
  {
    // get the node index vector from the mesh class:
    int ** nodeIndex = meshContainerPtr->getNodeIndexVector ();

    int ixMax = numMeshPointsX;
    int iyMax = numMeshPointsY;

    FILE *fp1 = fopen(filename,"w");

    for (j=0;j<iyMax;++j)
    {
      for (i=0;i<ixMax;++i)
      {
        int index = nodeIndex[i][j];
        fprintf(fp1,"%12.4e",xVec[index]  * scalingVars.x0);
        fprintf(fp1,"%12.4e",yVec[index]  * scalingVars.x0);
        fprintf(fp1,"%12.4e",VVec[index]  * scalingVars.V0);
        fprintf(fp1,"%12.4e",nnVec[index] * scalingVars.C0);
        fprintf(fp1,"%12.4e",npVec[index] * scalingVars.C0);
        fprintf(fp1,"%12.4e",CVec[index]  * scalingVars.C0);
        fprintf(fp1,"%12.4e",fabs(CVec[index]  * scalingVars.C0));
        fprintf(fp1,"%s","\n");
      }
      fprintf(fp1,"%s","\n");
    }

    fclose (fp1);
  }


  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::outputTxtData
// Purpose       : This function outputs data about each time step, that is
//                 easier to read as a text data, than as plotted data.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/15/03
//-----------------------------------------------------------------------------
bool Instance::outputTxtData ()
{
  bool bsuccess = true;

  int i;
  char filename[32];   for(i=0;i<32;++i) filename[i] = static_cast<char>(0);

  sprintf(filename,"%s_data.txt",outputName.c_str());

  double time = getSolverState().currTime_;

  FILE *fp1;
  if (callsOTXT<=0)
    fp1 = fopen(filename,"w");
  else
    fp1 = fopen(filename,"a");
  ++callsOTXT;

  // get maxs and mins.
  double VminOut   = 1.0e+99;
  double VmaxOut   =-1.0e+99;
  double NnMinOut  = 1.0e+99;
  double NnMaxOut  =-1.0e+99;
  double NpMinOut  = 1.0e+99;
  double NpMaxOut  =-1.0e+99;

  for (i=0;i<numMeshPoints; ++i)
  {
    if (VminOut  > ( VVec[i]*scalingVars.V0)) VminOut  = ( VVec[i]*scalingVars.V0);
    if (VmaxOut  < ( VVec[i]*scalingVars.V0)) VmaxOut  = ( VVec[i]*scalingVars.V0);
    if (NnMinOut > (nnVec[i]*scalingVars.C0)) NnMinOut = (nnVec[i]*scalingVars.C0);
    if (NnMaxOut < (nnVec[i]*scalingVars.C0)) NnMaxOut = (nnVec[i]*scalingVars.C0);
    if (NpMinOut > (npVec[i]*scalingVars.C0)) NpMinOut = (npVec[i]*scalingVars.C0);
    if (NpMaxOut < (npVec[i]*scalingVars.C0)) NpMaxOut = (npVec[i]*scalingVars.C0);
  }

  fprintf(fp1,"%s","\n");
  fprintf(fp1,"%s","---------------------------------------------------------\n");
  if (getSolverState().dcopFlag)
  {
    fprintf(fp1,"Global data for DC step %4d:\n", callsOTXT);
  }
  else
  {
    fprintf(fp1,"Global data for time step %4d:\n", callsOTXT);
  }
  fprintf(fp1,"Current Time = %12.4e\n" , time);
  fprintf(fp1,"       Vmin  = %12.4e\n", VminOut);
  fprintf(fp1,"       Vmax  = %12.4e\n", VmaxOut);
  fprintf(fp1,"       NnMin = %12.4e\n", NnMinOut);
  fprintf(fp1,"       NnMax = %12.4e\n", NnMaxOut);
  fprintf(fp1,"       NpMin = %12.4e\n", NpMinOut);
  fprintf(fp1,"       NpMax = %12.4e\n", NpMaxOut);

  fprintf(fp1,"%s","\n");

  // loop over the device interface nodes, sum the currents going into each one.
  std::vector<DeviceInterfaceNode>::iterator firstDI = dIVec.begin ();
  std::vector<DeviceInterfaceNode>::iterator lastDI  = dIVec.end   ();
  std::vector<DeviceInterfaceNode>::iterator iterDI;

  for (iterDI=firstDI; iterDI!=lastDI; ++iterDI)
  {
    // loop over the nodes of this device interface node:

     if ( !( meshContainerPtr->labelEdgeType (iterDI->eName) ) ) continue;

     mLabel * labelPtr = meshContainerPtr->getLabel(iterDI->eName);
     int firstNodeOfLabel = *(labelPtr->mNodeVector.begin());

     fprintf(fp1, "Information for electrode: %s\n", iterDI->eName.c_str());

     fprintf(fp1, "potential: %12.4e\n", scalingVars.V0* VVec[firstNodeOfLabel] );
     fprintf(fp1, "  current: %12.4e\n", iterDI->currentSum);
     fprintf(fp1, "  charge:  %12.4e\n", iterDI->chargeSum);
     fprintf(fp1, "  dIdVckt: %12.4e\n", iterDI->dIdVckt);
     fprintf(fp1, "  dQdVckt: %12.4e\n", iterDI->dQdVckt);
     fprintf(fp1, "%s","\n");
  }

  if (!calcConductanceFlag)
  {
    fprintf(fp1, "%s","NOTE:  The two-level Newton algorithm was not used.\n");
    fprintf(fp1, "%s","       This means that the conductances and capacitances\n");
    fprintf(fp1, "%s","       were not calculated.\n\n");
  }

  int iE1,iE2;

  fprintf(fp1, "%s","Conductance array: \n");
  fprintf(fp1,"%s","              ");
  for (iE2 = 0; iE2 < numElectrodes; ++iE2)
  {
    fprintf(fp1,"\t%14s",dIVec[iE2].eName.c_str());
  }
  fprintf(fp1,"%s","\n");

  for (iE1 = 0; iE1 < numElectrodes; ++iE1)
  {
    fprintf(fp1,"%14s",dIVec[iE1].eName.c_str());
    for (iE2 = 0; iE2 < numElectrodes; ++iE2)
    {
      fprintf(fp1,"\t%14.4e",condVec[iE1][iE2]);
    }
    fprintf(fp1,"%s","\n");
  }
  fprintf(fp1,"%s","\n");

  fprintf(fp1, "%s","Capacitance array: \n");
  fprintf(fp1,"%s","              ");
  for (iE2 = 0; iE2 < numElectrodes; ++iE2)
  {
    fprintf(fp1,"\t%14s",dIVec[iE2].eName.c_str());
  }
  fprintf(fp1,"%s","\n");

  for (iE1 = 0; iE1 < numElectrodes; ++iE1)
  {
    fprintf(fp1,"%14s",dIVec[iE1].eName.c_str());
    for (iE2 = 0; iE2 < numElectrodes; ++ iE2)
    {
      fprintf(fp1,"\t%14.4e",capVec[iE1][iE2]);
    }
    fprintf(fp1,"%s","\n");
  }
  fprintf(fp1,"%s","\n");

  fclose(fp1);

  return bsuccess;
}


//-----------------------------------------------------------------------------
// Function      : Instance::outputSgplot
// Purpose       : This function outputs a file that can be read by the
//                 plotting program sgplot.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 09/30/02
//-----------------------------------------------------------------------------
bool Instance::outputSgplot ()
{
  int i;

  // note: number of elements = N*Narray + Nvars.
  //   N  = number of nodes
  //   Narray = number of arrays
  //   Nvars  = number of non-array variables
  //   Nconst = number of constants

  UINT Narray = 6; // x,y,V,nn,np,C
  UINT Nconst = 3;

  RESHEAD resheadxyce =
  {   "SGFramework Result File  Version 1.0\n\x1a",  // text logo
      "@~!_RES" "ULT_!~@",          // signature
      "",                           // mesh file name
      Nconst,                       // number of constants
      0,                            // number of variables
      Narray,                       // number of arrays
      Narray*numMeshPoints,         // number of elements
      Narray,                       // number of 1D arrays
      0,                            // number of 2D arrays
      0,                            // number of 3D arrays
      1};                           // number of data sets

  RESHEAD *presheadxyce = &resheadxyce;

  // fix this:
  if (usingInternalMesh) meshFileName = outputName + ".msh";

  strcpy(presheadxyce->szMeshFile, meshFileName.c_str());

  // set up empty axlatconstxyce:
  XLATCONST axlatconstxyce[] =
  { { "NODES",     TYPE_ICONST, { numMeshPoints } },
    { "EDGES",     TYPE_ICONST, { numMeshEdges  } },
    { "TRIANGLES", TYPE_ICONST, { numMeshCells  } }
  };

  DAXLATARRAY AxlatArray(numMeshPoints);  // see the constructor

  XLATARRAY *paxlat = AxlatArray.GetPointer(0);

  char filename[32];   for(i=0;i<32;++i) filename[i] = static_cast<char>(0);
  sprintf(filename,"%s_%03d.res",outputName.c_str(),callsOSG);
  ++callsOSG;

  FILE *nHandle = fopen(filename, "w");

  UINT cConst = Nconst;
  // XLATARRAY union fix makes this unneeded
  // for(i = 0; i < Nconst; ++i)
  //   axlatconstxyce[i].data.n = static_cast<int> (axlatconstxyce[i].data.r);

  fwrite(presheadxyce, sizeof(RESHEAD),1,nHandle);
  if (cConst) fwrite(axlatconstxyce, sizeof(XLATCONST), cConst, nHandle);
  if (Narray ) fwrite(        paxlat, sizeof(XLATARRAY), Narray, nHandle);

  // In Kevin's original setup, these two lines are in place.
  // The reason for this is so that he can keep adding data sets
  // to the same file.  To do so, he has to update the cSset
  // attribute of RESHEAD, but nothing else.  So he rewrites that,
  // then goes back to the end of the file.
  // Anyway, it isn't neccessary here.  I only put one data set in
  // any file.
  //presheadxyce->cSet++;
  //fwrite(presheadxyce, sizeof(RESHEAD),1,nHandle);
  //fseek(nHandle, 0L, SEEK_END);

  // output all the data arrays:
  for (i=0;i<numMeshPoints;++i) outputVec[i] = scalingVars.x0 * xVec[i];
  fwrite(  &(outputVec[0]), sizeof(double), numMeshPoints, nHandle);

  for (i=0;i<numMeshPoints;++i) outputVec[i] = scalingVars.x0 * yVec[i];
  fwrite(  &(outputVec[0]), sizeof(double), numMeshPoints, nHandle);

  for (i=0;i<numMeshPoints;++i) outputVec[i] = scalingVars.V0 * VVec[i];
  fwrite(  &(outputVec[0]), sizeof(double), numMeshPoints, nHandle);

  for (i=0;i<numMeshPoints;++i) outputVec[i] = scalingVars.C0 * nnVec[i];
  fwrite( &(outputVec[0]), sizeof(double), numMeshPoints, nHandle);

  for (i=0;i<numMeshPoints;++i) outputVec[i] = scalingVars.C0 * npVec[i];
  fwrite( &(outputVec[0]), sizeof(double), numMeshPoints, nHandle);

  for (i=0;i<numMeshPoints;++i) outputVec[i] = scalingVars.C0 * CVec[i];
  fwrite(  &(outputVec[0]), sizeof(double), numMeshPoints, nHandle);

  fclose(nHandle);

  return true;
}

} // namespace TwoDPDE
} // namespace Device
} // namespace Xyce

//-----------------------------------------------------------------------------
// Function      : DAXLATARRAY::set
// Purpose       : This function (and class) is only used in the context
//                 of outputting sgplot files.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 09/30/02
//-----------------------------------------------------------------------------
void DAXLATARRAY::set (const char *name, UINT uOffset, UINT cDim, UINT ac0, UINT ac1, UINT ac2)
{
  XLATARRAY xlattmp;
  char tmpname[LEN_IDENT+1];
  int istrl = strlen(name);
  int imax;

  if(istrl >= LEN_IDENT-1)  imax = LEN_IDENT-1;
  else        imax = istrl;

  int i;
  for(i=0;i<LEN_IDENT+1;++i) tmpname[i] = 0;

  for(i=0;i<imax;++i)
    tmpname[i] = name[i];

  sprintf(xlattmp.szName,"%s",tmpname);

  xlattmp.uOffset = uOffset;
  xlattmp.cDim    = cDim;
  xlattmp.acElements[0] = ac0;
  xlattmp.acElements[1] = ac1;
  xlattmp.acElements[2] = ac2;
  Add(xlattmp);

}

//-----------------------------------------------------------------------------
// Function      : DAXLATARRAY::DAXLATARRAY
// Purpose       : constructor
// Special Notes : This function (and class) is only used in the context
//                 of outputting sgplot files.
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 09/30/02
//-----------------------------------------------------------------------------
DAXLATARRAY::DAXLATARRAY(int numMeshPoints)
{
  int c = 10;
  cElements = 0;
  uInc = cSize = c;
  dynarray = new XLATARRAY[c];

  char tmpArg[16]; for(int i=0;i<16;++i) tmpArg[i] = 0;

  // set up axlatarrayicp:
  set( "x" , 0*numMeshPoints,1,numMeshPoints,0,0);
  set( "y" , 1*numMeshPoints,1,numMeshPoints,0,0);
  set( "V" , 2*numMeshPoints,1,numMeshPoints,0,0);
  set( "Ne", 3*numMeshPoints,1,numMeshPoints,0,0);
  set( "Np" , 4*numMeshPoints,1,numMeshPoints,0,0);
  set( "C" , 5*numMeshPoints,1,numMeshPoints,0,0);
}


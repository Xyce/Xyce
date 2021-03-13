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

//-----------------------------------------------------------------------------
//
// Purpose        : Specification of Warped MPDE phase condition (graph & load)
//
// Special Notes  :
//
// Creator        : Todd Coffey, SNL, 1414
//
// Creation Date  : 12/7/06
//
//
//
//
//-----------------------------------------------------------------------------

#include <Xyce_config.h>


#include <N_MPDE_WarpedPhaseCondition.h>
#include <N_PDS_Comm.h>

using std::max;

Teuchos::RCP<std::vector<int> > N_MPDE_WarpedPhaseCondition::getPhaseGraph() const
{
  if ((warpMPDEOSCOUT_ == -1) && (warpPhase_ != 0))
  {
    Xyce::Report::DevelFatal0().in("N_MPDE_WarpedPhaseCondition::getPhaseGraph")
      << " No value for oscout which is required by specified phase equation";
  }
  Teuchos::RCP<std::vector<int> > phaseGraph = Teuchos::rcp(new std::vector<int>);
  if ( warpPhase_ == 0 )
  {
    // Add graph for phase condition equation for MPDE:  omega-1=0
    phaseGraph->push_back(omegaGID_);
  }
  else if ( warpPhase_ == 1 )
  {
    // Add graph for phase condition 1 for WaMPDE:\hat{x}_w(t_1,0) - alpha = 0
    phaseGraph->push_back(warpMPDEOSCOUT_);
  }
  else if ( warpPhase_ == 2 )
  {
    // Add graph for phase condition 2 for WaMPDE:
    // \omega*((\hat{x}_w(t_1,0)-hat{x}_w(t_1,-h2))/h2 - alpha) = 0
    phaseGraph->push_back( warpMPDEOSCOUT_ );
    phaseGraph->push_back( warpMPDEOSCOUT_ + (size_-1)*offset_ );
    phaseGraph->push_back( omegaGID_ );
  }
  else if ( warpPhase_ == 3 )
  {
    // Add graph for phase condition 2 for WaMPDE:
    // \omega*((\hat{x}_w(t_1,h2)-hat{x}_w(t_1,-h2))/(2*h2) - alpha) = 0
    phaseGraph->push_back( warpMPDEOSCOUT_ + (   1   )*offset_ );
    phaseGraph->push_back( warpMPDEOSCOUT_ + (size_-1)*offset_ );
    phaseGraph->push_back( omegaGID_ );
  }
  else
  {
    Xyce::Report::UserWarning()
      << " Unrecognized value for WaMPDE Phase option";
  }
  return(phaseGraph);
}


double
N_MPDE_WarpedPhaseCondition::getPhaseCondition(
  Xyce::Linear::BlockVector *   bXptr,
  std::vector<double> & fastTimes) const
{
  Xyce::Linear::BlockVector& bX = *bXptr;
  int BlockCount = bX.blockCount();
  double phaseValue = 0.0, tmpPhaseValue = 0.0;
  int omegaLID = bX.pmap()->globalToLocalIndex(omegaGID_);
 
  // Get omega to all processors since it is not certain that the omega node is
  // on the same processor as the warpMPDEOSCOUT node.
  double omega = 0.0, tmpOmega = 0.0;
  if (omegaLID >= 0)
  {
    tmpOmega = bX[omegaLID];
  }
  bX.pmap()->pdsComm().sumAll( &tmpOmega, &omega, 1 );

  // Get the local ID for the warpMPDEOSCOUT node, assume all block shifts are also
  // on the same processor.
  int warpMPDEOSCOUTLID = bX.pmap()->globalToLocalIndex(warpMPDEOSCOUT_);
  if (warpMPDEOSCOUTLID >= 0)
  {
    if ( warpPhase_ == 0 )
    {
      tmpPhaseValue = omega - 1.0;
    }
    else if ( warpPhase_ == 1 )
    {
      // Add graph for phase condition 1 for WaMPDE:  \hat{x}_w(t_1,0) - alpha = 0
      tmpPhaseValue = bX[warpMPDEOSCOUTLID] - warpPhaseCoeff_;
    }
    else if ( warpPhase_ == 2 )
    {
      // Add graph for phase condition 2 for WaMPDE:
      // \omega*((\hat{x}_w(t_1,0)-hat{x}_w(t_1,-h2))/h2 - alpha) = 0
      int shift = (BlockCount-1)*offset_;
      int xwmLID = bX.pmap()->globalToLocalIndex(warpMPDEOSCOUT_+shift);
      double xwp = bX[warpMPDEOSCOUTLID];
      double xwm = bX[xwmLID];
      double invh2 = 1.0 / (fastTimes[BlockCount] - fastTimes[BlockCount-1]);
      tmpPhaseValue = omega*(xwp - xwm)*invh2 - warpPhaseCoeff_ ;
    }
    else if ( warpPhase_ == 3 )
    {
      // Add graph for phase condition 2 for WaMPDE:
      // \omega*((\hat{x}_w(t_1,h2)-hat{x}_w(t_1,-h2))/(2*h2) - alpha) = 0
      int shiftp = offset_;
    int shiftm = (BlockCount-1)*offset_;
    int xwpLID = bX.pmap()->globalToLocalIndex(warpMPDEOSCOUT_+shiftp);
    int xwmLID = bX.pmap()->globalToLocalIndex(warpMPDEOSCOUT_+shiftm);
    double xwp = bX[xwpLID];
    double xwm = bX[xwmLID];
    double max_xw = max(abs(xwp),abs(xwm));
    double invh2 = 1.0 / ((fastTimes[BlockCount] - fastTimes[BlockCount-1]) + (fastTimes[1] - fastTimes[0]));
      tmpPhaseValue = omega*(xwp - xwm)*(1.0 / max_xw) - warpPhaseCoeff_ ;
    }
    else
    {
      Xyce::Report::UserWarning()
        << " Unrecognized value for WaMPDE Phase option";
    }
  }
 
  // Communicate to all processors.
  bX.pmap()->pdsComm().sumAll( &tmpPhaseValue, &phaseValue, 1 );
   
  return(phaseValue);
}

void
N_MPDE_WarpedPhaseCondition::getPhaseConditionDerivative(
  Xyce::Linear::BlockVector *   bXptr,
  std::vector<double> & fastTimes,
  std::vector<int> *    colIndicesPtr,
  std::vector<double> * coeffsPtr) const
{
  Xyce::Linear::BlockVector & bX = *bXptr;
  std::vector<int> & colIndices = *colIndicesPtr;
  std::vector<double> & coeffs = *coeffsPtr;
  int BlockCount = bX.blockCount();
  
  // Get omega to all processors since it is not certain that the omega node is
  // on the same processor as the warpMPDEOSCOUT node.
  double omega = 0.0, tmpOmega = 0.0;
  int omegaLID = bX.pmap()->globalToLocalIndex(omegaGID_);
  if (omegaLID >= 0)
  {
    tmpOmega = bX[omegaLID];
  }
  bX.pmap()->pdsComm().sumAll( &tmpOmega, &omega, 1 );

  int warpMPDEOSCOUTLID = bX.pmap()->globalToLocalIndex(warpMPDEOSCOUT_);
  if ( warpMPDEOSCOUTLID >= 0 )
  {
    if ( warpPhase_ == 0 )
    {
      // put derivative wrt x_v of (omega-1) = 0 into omegaGID_,x_v locations
      coeffs.resize(1);
      coeffs[0] = 1.0;
      colIndices.resize(1);
      colIndices[0] = omegaGID_;
    }
    else if ( warpPhase_ == 1 )
    {
      // put derivative wrt x_v of (x_v-warpMPDECoeff) = 0 into omegaGID_,x_v location
      coeffs.resize(1);
      coeffs[0] = 1.0;
      colIndices.resize(1);
      colIndices[0] = warpMPDEOSCOUT_;
    }
    else if ( warpPhase_ == 2 )
    {
      // put derivative wrt x_v, xv+shift, and omega of
      // \omega*((\hat{x}_w(t_1,0)-\hat{x}_w(t_1,-h2))/h2 - alpha)
      //   = omega/h2 into  omegaGID_,x_v location
      //   = -omega/h2 into omegaGID_,x_v+shift location
      //   = (\hat{x}_w(t_1,0)-\hat{x}_w(t_1,-h2))/h2 into omegaGID_,omegaGID_ location
      int shift = (BlockCount-1)*offset_;
      int xwmLID = bX.pmap()->globalToLocalIndex(warpMPDEOSCOUT_+shift);
      double xwp = bX[warpMPDEOSCOUTLID];
      double xwm = bX[xwmLID];
      double invh2 = fastTimes[BlockCount] - fastTimes[BlockCount-1];
      coeffs.resize(3);
      coeffs[0] =  omega*invh2;
      coeffs[1] = -omega*invh2;
      coeffs[2] = (xwp - xwm)*invh2;
      colIndices.resize(3);
      colIndices[0] = warpMPDEOSCOUT_;
      colIndices[1] = warpMPDEOSCOUT_+shift;
      colIndices[2] = omegaGID_;
    }
    else if ( warpPhase_ == 3 )
    {
      // put derivative wrt x_v, xv+shift, and omega of
      // \omega*((\hat{x}_w(t_1,h2)-hat{x}_w(t_1,-h2))/(2*h2) - alpha)
      //   = omega/(2*h2) into  omegaGID_,x_v+shiftp location
      //   = -omega/(2*h2) into omegaGID_,x_v+shiftm location
      //   = (\hat{x}_w(t_1,h2)-\hat{x}_w(t_1,-h2))/(2*h2) into omegaGID_,omegaGID_ location
      int shiftp = offset_;
      int shiftm = (BlockCount-1)*offset_;
      int xwpLID = bX.pmap()->globalToLocalIndex(warpMPDEOSCOUT_+shiftp);
      int xwmLID = bX.pmap()->globalToLocalIndex(warpMPDEOSCOUT_+shiftm);
      double xwp = bX[xwpLID];
      double xwm = bX[xwmLID];
      double invh2 = (fastTimes[BlockCount] - fastTimes[BlockCount-1]) + (fastTimes[1] - fastTimes[0]);
      coeffs.resize(3);
      coeffs[0] = omega*0.5*invh2;
      coeffs[1] = -omega*0.5*invh2;
      coeffs[2] = (xwp-xwm)*0.5*invh2;
      colIndices.resize(3);
      colIndices[0] = warpMPDEOSCOUT_+shiftp;
      colIndices[1] = warpMPDEOSCOUT_+shiftm;
      colIndices[2] = omegaGID_;
    }
    else
    {
      Xyce::Report::UserWarning()
        << " Unrecognized value for WaMPDE Phase option\n";
    }
  }
  else 
  {
    coeffs.resize(0);
  }

#ifdef Xyce_PARALLEL_MPI
  // Now get this information to the rest of the processors.
  int numCoeffs = 0, tmpNumCoeffs = coeffs.size();
  bX.pmap()->pdsComm().maxAll( &tmpNumCoeffs, &numCoeffs, 1 );
  std::vector<int> tmpColIndices(numCoeffs, -1);
  std::vector<double> tmpCoeffs(numCoeffs, 0.0);
  if ( warpMPDEOSCOUTLID < 0 )
  {
    coeffs.resize(numCoeffs);
    colIndices.resize(numCoeffs);
  }
  else
  {
    tmpColIndices = colIndices;
    tmpCoeffs = coeffs;
  }
  bX.pmap()->pdsComm().sumAll( &tmpCoeffs[0], &coeffs[0], numCoeffs );
  bX.pmap()->pdsComm().maxAll( &tmpColIndices[0], &colIndices[0], numCoeffs );
#endif 
}

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

//-----------------------------------------------------------------------------
// Purpose       : 
// Special Notes :
// Creator       : Robert Hoekstra, 9233, Computational Sciences
// Creation Date : 5/6/04
//
//-----------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------

#include <iostream>

// ----------   Xyce Includes   ----------
#include <N_ERH_ErrorMgr.h>
#include <N_MPDE_Discretization.h>
#include <N_UTL_LogStream.h>
#include <N_UTL_FeatureTest.h>

using Xyce::DEBUG_MPDE;

//-----------------------------------------------------------------------------
// Function      : N_MPDE_Discretization::N_MPDE_Discretization
// Purpose       : Constructor
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, 9233, Computational Sciences
// Creation Date : 5/6/04
//-----------------------------------------------------------------------------
N_MPDE_Discretization::N_MPDE_Discretization( Type type, int order )
  : type_(type),
    order_(order)
{
  switch( type )
  {
    case Backward:
      switch( order )
      {
        case 1:
          start_ = -1;
          width_ =  2;
          coeffs_.resize( width_ );
          coeffs_[0] = -1.0;
          coeffs_[1] =  1.0;
          break;
          
        case 2:
          start_ = -2;
          width_ =  3;
          coeffs_.resize( width_ );
          coeffs_[0] =  1.0;
          coeffs_[1] = -4.0;
          coeffs_[2] =  3.0;
          break;
          
        case 3:
          start_ = -3;
          width_ =  4;
          coeffs_.resize( width_ );
          coeffs_[0] = -1.0;
          coeffs_[1] =  4.5;
          coeffs_[2] = -9.0;
          coeffs_[3] =  5.5;
          break;
          
        default:
          Xyce::Report::UserFatal()
            << "MPDE Discretization Error.  Backward differences only supported for order=1, 2 or 3.";
      }
      break;
      
    case Centered:
      switch( order )
      {
        case 1:
          Xyce::Report::UserWarning()
            << "MPDE Discretization Warning.  Central differences requested with order = 1.  Defaulting to order = 2.";
          order_ = 2;
            
        case 2:
          start_ = -1;
          width_ =  3;
          coeffs_.resize( width_ );
          coeffs_[0] = -1.0;
          coeffs_[1] =  0.0;
          coeffs_[2] =  1.0;
          break;
          
        case 3:
          start_ = -2;
          width_ =  5;
          coeffs_.resize( width_ );
          coeffs_[0] =  1.0/3.0;
          coeffs_[1] = -8.0/3.0;
          coeffs_[2] =  0.0;
          coeffs_[3] =  8.0/3.0;
          coeffs_[4] = -1.0/3.0;
          break;
          
        default:
          Xyce::Report::UserFatal()
            <<  "MPDE Discretization Error.  Central differences only supported for order=2 and 3.";
      
      }
      break;
    
    case Forward:
      switch( order )
      {
        case 1:
          start_ =  0;
          width_ =  2;
          coeffs_.resize( width_ );
          coeffs_[0] = -1.0;
          coeffs_[1] =  1.0;
          break;

        case 2:
          start_ =  0;
          width_ =  3;
          coeffs_.resize( width_ );
          coeffs_[0] = -3.0;
          coeffs_[1] =  4.0;
          coeffs_[2] = -1.0;
          break;
          
        case 3:
          start_ =  0;
          width_ =  4;
          coeffs_.resize( width_ );
          coeffs_[0] = -5.5;
          coeffs_[1] =  9.0;
          coeffs_[2] = -4.5;
          coeffs_[3] =  1.0;
          break; 
          
        default:
          Xyce::Report::UserFatal()
            << "MPDE Discretization Error.  Forward differences only supported for order=1, 2 and 3.";
      }
      break;
      
    default:
      Xyce::Report::UserFatal() 
        << "MPDE Discretization Error.  Unspecified differencing scheme for MPDE fast time scale.";
      break;
  }
  
  /* 
     GenerateCoeffs_ assumes even spacing between grid points.  This
     isn't a problem for first order methods, but for higher order 
     methods, the denominator which can be h1 + h2 becomes 2h and the
     prefactor, 2 here, is included with the coefficients.  This 
     causes problems when trying to apply the coefficients to a non-
     uniform grid as one must calculate what this prefactor was and 
     factor it out.
  
     So for now we'll just statically generate first, second and third
     order coefficients for now
  */
  
  // GenerateCoeffs_( type, order, coeffs_ );
  
  if (DEBUG_MPDE)
  {
    Xyce::dout() << "MPDE Fast Time Disc\n";
    Xyce::dout() << "-------------------\n";
    if( type == Backward ) Xyce::dout() << "Disc: Backward Difference\n";
    else if ( type == Centered ) Xyce::dout() << "Disc: Centered Difference\n";
    else if ( type == Forward ) Xyce::dout() << "Disc: Forward Difference\n";

    Xyce::dout() << "Order: " << order_ << std::endl;
    Xyce::dout() << "Start: " << start_ << std::endl;
    Xyce::dout() << "Width: " << width_ << std::endl;
    Xyce::dout() << "Coefficients: ";
    for( int i = 0; i < width_; ++i ) Xyce::dout() << coeffs_[i] << " ";
    Xyce::dout() << std::endl;
    Xyce::dout() << "-------------------\n";
  }
}

//-----------------------------------------------------------------------------
// Function      : N_MPDE_Discretization::GenerateCoeffs_
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, 9233, Computational Sciences
// Creation Date : 5/6/04
//-----------------------------------------------------------------------------
void N_MPDE_Discretization::GenerateCoeffs_( Type type, int order, std::vector<double> & coeffs )
{
  /*
   * Function to calculate finite difference weights for derivatives of
   * order 0 through m using n+1 grid points.
   *
   * diff_flag          the type of stencil to use, either FORWARD_DIFFERENCE,
   *                    BACKWARD_DIFFERENCE, or CENTERED_DIFFERENCE.
   * m                  the largerst degree of derivative for which weights
   *                    are sought.
   * n                  one less than the number of grid points (related to
   *                    the order of accuracy of the difference stencil)
   * c                  a m+1 by n+1 array of finite difference weights.
   *                    c[m][0..n] are the weights for the mth derivative.
   *
   * This algorithm is taken almost verbatim by Eric Phipps from:
   * B. Fornberg, Calculation of Finite Difference Formulas, SIAM Rev. 
   * Vol. 40, No. 3, pp 685-691, Sept. 1998.
   */
  double c1=1.0,c2=0.0,c3=0.0,c4=0.0,c5=0.0;
  int i,j,k,mn;

  int m = 1;
  int n = order;

  switch (type)
  {
   case Backward:
    c4 = static_cast<double> (-n);
    break;
   case Centered:
    c4 = -0.5* static_cast<double> (n);
    break;
   case Forward:
    // do nothing
    break;
  }
  
  std::vector< std::vector<double> > c(2);
  c[0].resize(n+1);
  c[1].resize(n+1);
  for (k=0; k<=m; ++k)
    for (j=0; j<=n; ++j)
      c[k][j] = 0.;
  c[0][0] = 1.;
  
  for (i=1; i<=n; ++i)
  {
    mn = i<=m ? i : m;
    c2 = 1.0;
    c5 = c4;
    
    switch(type)
    {
     case Backward:
      c4 = static_cast<double> (-n+i);
      break;
     case Centered:
      c4 = -0.5 * static_cast<double> (n) + static_cast<double> (i);
      break;
     case Forward:
      // do nothing
      break;
    }
    
    for (j=0; j<=i-1; ++j)
    {
      c3 = static_cast<double> (i-j);
      c2 = c2*c3;
      if (j == i-1) {
        for (k = mn; k>=1; k--)
          c[k][i] = c1*(k*c[k-1][i-1]-c5*c[k][i-1])/c2;
        c[0][i] = -c1*c5*c[0][i-1]/c2;
      }
      for (k=mn; k>=1; k--)
        c[k][j] = (c4*c[k][j] - static_cast<double> (k) * c[k-1][j])/c3;
      c[0][j] = c4*c[0][j]/c3;
    }
    c1 = c2;
  }
      
  coeffs.resize(n+1);
  for( i = 0; i < n+1; ++i ) coeffs[i] = c[1][i];
}


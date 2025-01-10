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

//-----------------------------------------------------------------------------
//
// Purpose        : Interpolator classes
//
// Special Notes  :
//
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 01/31/12
//-----------------------------------------------------------------------------

#ifndef Xyce_N_UTL_Interpolators_h
#define Xyce_N_UTL_Interpolators_h

#include <Sacado_No_Kokkos.hpp>

#include <complex>

#include <N_ERH_ErrorMgr.h>

namespace Xyce {
namespace Util {

template <typename ScalarA, typename ScalarB>
inline bool greaterThan(ScalarA & left, ScalarB & right) { return (left > right); }

template <>
inline bool greaterThan(std::complex<double> & left, std::complex<double> & right) { return (std::real(left) > std::real(right)); }

template <>
inline bool greaterThan(const std::complex<double> & left, const std::complex<double> & right) { return (std::real(left) > std::real(right)); }


//-----------------------------------------------------------------------------
// Class         : interpolator base class
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter, SNL
// Creation Date : 1/31/2012
//-----------------------------------------------------------------------------
template <typename ScalarT>
class interpolator
{
public:
  interpolator (){};
  virtual ~interpolator (){};

  virtual void clear (){};

  virtual void init (const std::vector<ScalarT> & xa, const std::vector<ScalarT> & ya){};

  virtual void eval (
     const std::vector<ScalarT> & xa, const std::vector<ScalarT> & ya,
     const ScalarT & x, ScalarT & y) const {};

  virtual void evalDeriv (
     const std::vector<ScalarT> & xa, const std::vector<ScalarT> & ya,
     const ScalarT & x, ScalarT & dydx) const {};

  virtual void evalDeriv2 (
     const std::vector<ScalarT> & xa, const std::vector<ScalarT> & ya,
     const ScalarT & x, ScalarT & ypp) const {};

  virtual void evalInteg (
     const std::vector<ScalarT> & xa, const std::vector<ScalarT> & ya,
     const ScalarT & a, const ScalarT & b, ScalarT & result) const {};

  virtual void
    getCoefs (
     const std::vector<ScalarT> & xa, const std::vector<ScalarT> & ya,
        size_t index, std::vector<ScalarT> & coefs) {};

  inline size_t
  binarySearch(
     const std::vector<ScalarT> & xa,
     const ScalarT & x,
     size_t index_lo,
     size_t index_hi) const;

  inline ScalarT
  integ_eval (
     const ScalarT & ai, const ScalarT & bi, const ScalarT & ci,
     const ScalarT & di, const ScalarT & xi, const ScalarT & a,
     const ScalarT & b) const;
};

//-----------------------------------------------------------------------------
// Function      : interpolator<ScalarT>::binarySearch
// Purpose       : Perform a binary search of an array of values.
// Special Notes : 
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 01/31/2012
// ----------------------------------------------------------------------------
template <typename ScalarT>
inline size_t
interpolator<ScalarT>::binarySearch(
   const std::vector<ScalarT> & xa,
   const ScalarT & x,
   size_t index_lo,
   size_t index_hi) const
{
  size_t ilo = index_lo;
  size_t ihi = index_hi;
  while(ihi > ilo + 1)
  {
    //size_t i = (ihi + ilo)/2;
    size_t i = (ihi + ilo) >> 1;

    if(greaterThan(xa[i],x))
    {
      ihi = i;
    }
    else
    {
      ilo = i;
    }
  }
  return ilo;
}

//-----------------------------------------------------------------------------
// Function      : interpolator<ScalarT>::integ_eval
//
// Purpose       : function for doing the spline integral evaluation
//                 which is common to most of the interpolators.
//
// Special Notes : 
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 01/31/2012
// ----------------------------------------------------------------------------
template <typename ScalarT>
inline ScalarT
interpolator<ScalarT>::integ_eval (
   const ScalarT & ai,
   const ScalarT & bi,
   const ScalarT & ci,
   const ScalarT & di,
   const ScalarT & xi,
   const ScalarT & a,
   const ScalarT & b) const
{
  const ScalarT r1 = a - xi;
  const ScalarT r2 = b - xi;
  const ScalarT r12 = r1 + r2;
  const ScalarT bterm = 0.5 * bi * r12;
  const ScalarT cterm = (1.0 / 3.0) * ci * (r1 * r1 + r2 * r2 + r1 * r2);
  const ScalarT dterm = 0.25 * di * r12 * (r1 * r1 + r2 * r2);
  return (b - a) * (ai + bterm + cterm + dterm);
}

//-----------------------------------------------------------------------------
// Class         : akima spline class
// Purpose       : 
// Special Notes : References:
//
// The original reference for this type of spline is in this paper:
//
//   Hiroshi Akima. 1970. A New Method of Interpolation and Smooth Curve Fitting 
//   Based on Local Procedures. J. ACM 17, 4 (October 1970), 589-602. 
//   DOI=http://dx.doi.org/10.1145/321607.321609
//
// There are several variants of the Akima spline found in the literature.   
// However, the version used in this class matches that of the original 1970 paper.
//
// A different variation of the akima spline, which has rounded corners, is 
// implemented in the "wodicka" class.  See the comments in that class for details.
//
// Creator       : Eric Keiter, SNL
// Creation Date : 1/31/2012
//-----------------------------------------------------------------------------
template <typename ScalarT>
class akima: public interpolator<ScalarT>
{
public:
  akima () {};
  ~akima () { clear(); };

  void init ( const std::vector<ScalarT> & xa, const std::vector<ScalarT> & ya);

  void clear () { p1.clear(); p2.clear(); p3.clear(); m.clear(); };

  void eval (
     const std::vector<ScalarT> & xa, const std::vector<ScalarT> & ya,
     const ScalarT & x, ScalarT & y) const;

  void evalDeriv (
     const std::vector<ScalarT> & xa, const std::vector<ScalarT> & ya,
     const ScalarT & x, ScalarT & dydx) const;

  void evalDeriv2 (
     const std::vector<ScalarT> & xa, const std::vector<ScalarT> & ya,
     const ScalarT & x, ScalarT & ypp) const;

  void evalInteg (
     const std::vector<ScalarT> & xa, const std::vector<ScalarT> & ya,
     const ScalarT & a, const ScalarT & b, ScalarT & result) const;

  void getCoefs (
     const std::vector<ScalarT> & xa, const std::vector<ScalarT> & ya,
        size_t index, std::vector<ScalarT> & coefs)
    {
      coefs.clear();
      if (index < xa.size())
      {
        coefs.resize(4,0.0);
        coefs[0] = ya[index];
        coefs[1] = p1[index];
        coefs[2] = p2[index];
        coefs[3] = p3[index];
      }
    };

public:
  std::vector<ScalarT>  p1;
  std::vector<ScalarT>  p2;
  std::vector<ScalarT>  p3;
  std::vector<ScalarT>  m;
};

//-----------------------------------------------------------------------------
// Function      : akima<ScalarT>::init
// Purpose       :
// Special Notes : 
//
// Note that the 'm' indices are offset by 2 on either side of the array,
// compared with the x indices.  ie, m[2] corresponds to x[0], etc.  This is 
// because the Akima algorithm requires 5 points for a given data 
// point - two below and two above.  At each boundary, two additional 
// artificial points which approximate the edge-based derivatives had to be 
// added.
//
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 1/31/2012
// ----------------------------------------------------------------------------
template <typename ScalarT>
void akima<ScalarT>::init (
   const std::vector<ScalarT> & xa,
   const std::vector<ScalarT> & ya)
{
  size_t size = xa.size();
  size_t i;

  if (size<=0)
  {
    Report::DevelFatal().in("akima<ScalarT>::init") << "Array size  = " << size << ".  Inteprolation failed";
  }

  if (p1.size() != size) p1.resize(size,0.0);
  if (p2.size() != size) p2.resize(size,0.0);
  if (p3.size() != size) p3.resize(size,0.0);

  // m represents edges between points.  So, original size of m is n-1, if the 
  // size of xa and ya is n.  Therefore, when adding the additional 4 points to 
  // m, required by the method, the size of m should be n+3.
  if (m.size() != size+3) 
  {
    m.resize(size+3,0.0);
  }

  // set up the edge derivatives
  if (xa.size() > 1)
  {
    for (i = 0; i < size-1; i++)
    {
      m[i+2] = (ya[i + 1] - ya[i]) / (xa[i + 1] - xa[i]);
    }
  }

  // As described in the original paper, boundary conditions are handled slightly differently.  
  // The slopes for two extra segments on either side of the array need to be estimated.
  //
  // For a non-periodic curve:
  //
  //    indices from 1970 paper:                  indices used here:
  //   m_{-1}  = 3*m_{1} − 2*m_{2}      -> m_{0}   = 3*m_{2} − 2*m_{3}   (adding one)
  //   m_{0}   = 2*m_{1} − m_{2}        -> m_{1}   = 2*m_{1} − m_{2}     (adding one)
  //   m_{n}   = 2*m_{n−1} − m_{n−2}    -> m_{n+1} = 2*m_{n} − m_{n−1}   (adding one)
  //   m_{n+1} = 3*m_{n−1} − 2*m_{n−2}  -> m_{n+2} = 3*m_{n] - 2*m_{n-1} (adding one)
  //
  // For a periodic curve: 
  //
  //    indices from 1970 paper:          indices used here:
  //   m_{−1}  = m_{n−2}        ->  m_0 = m_{n-1}
  //   m_{0}   = m_{n−1}        ->  m_1 = m_{n}
  //   m_{n}   = m_{1}          ->  m_{n+1} = m_{2}
  //   m_{n+1} = m_{2}          ->  m_{n+2} = m_{3}    
  //
  // Assuming non-periodic boundary conditions:
  m[0] = 3.0 * m[2] - 2.0 * m[3];
  m[1] = 2.0 * m[2] - m[3];
  m[size + 1] = 2.0 * m[size] - m[size-1];
  m[size + 2] = 3.0 * m[size] - 2.0 * m[size-1];

  // Periodic boundary conditions:
  //m[0] = m[size-1]
  //m[1] = m[size]
  //m[size + 1] =  m[2]
  //m[size + 2] =  m[3]
 
  ScalarT t1 = 0.0;
  ScalarT t2 = 0.0;

  // Original Akima 1970 paper:
  std::vector<ScalarT> t(size,0.0);
  for (i = 0; i < size; i++)
  {
    const ScalarT dm32 = std::abs(m[i + 3] - m[i + 2]);
    const ScalarT dm10 = std::abs(m[i + 1] - m[i]);

    if ((dm32+dm10) == 0.0)
    {
      t[i] = 0.5*(m[i+1] + m[i+2]);
    }
    else
    {
      t[i] = (dm32 * m[i+1] + dm10 * m[i+2])/(dm32 + dm10);
    }
  }

  for (i = 0; i < size-1; i++)
  {
    t1 = t[i];
    t2 = t[i+1];

    // see equations 3-6 of the original Akima 1970 paper.  
    // (p0 is used implicitly in the eval function, so not computed/stored)
    const ScalarT dx = xa[i + 1] - xa[i];
    p1[i] = t1;
    p2[i] = (3.0 * m[i+2] - 2.0 * t1 - t2) / dx;
    p3[i] = (t1 + t2 - 2.0 * m[i+2]) / (dx * dx);
  }
}

//-----------------------------------------------------------------------------
// Function      : akima<ScalarT>::eval
// Purpose       :
// Special Notes : 
//
//  y(x) = p0+ pl(x-xl) + p2(x-x1)^2+ p3(x-x1)^3
//
//  where:
//
//  p0 = y1 = y(x1)  (result of the binary search)
//  p1 = t1 - estimated slope dydx at x1, based on equation 1 from original paper
//  p2 = [3(y2 - y,)/(x2 - xl) - 2tl - t2]/(x2 - xl)
//  p3 = [t1 + t2 - 2(y2 - yl)/(x2 - xl)]/(x2 - xl)^2
//
//  p0,p1,p2 and p3 were all computed in the "init" function.
//
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 1/31/2012
// ----------------------------------------------------------------------------
template <typename ScalarT>
void akima<ScalarT>::eval (
   const std::vector<ScalarT> & xa,
   const std::vector<ScalarT> & ya,
   const ScalarT & x,
   ScalarT & y) const
{
  size_t size = xa.size();
  size_t index = this->binarySearch (xa, x, 0, size - 1);

  const ScalarT delx = x - xa[index];
  y = ya[index] + delx * (p1[index] + delx * (p2[index] + p3[index] * delx));
  return;
}

//-----------------------------------------------------------------------------
// Function      : akima<ScalarT>::evalDeriv
// Purpose       :
// Special Notes : 
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 1/31/2012
// ----------------------------------------------------------------------------
template <typename ScalarT>
void akima<ScalarT>::evalDeriv (
   const std::vector<ScalarT> & xa,
   const std::vector<ScalarT> & ya,
   const ScalarT & x,
   ScalarT & dydx) const
{
  size_t size = xa.size();
  size_t index = this->binarySearch (xa, x, 0, size - 1);

  ScalarT delx = x - xa[index];
  dydx = p1[index] + delx * (2.0 * p2[index] + 3.0 * p3[index] * delx);
  return;
}

//-----------------------------------------------------------------------------
// Function      : akima<ScalarT>::evalDeriv2
// Purpose       :
// Special Notes : 
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 1/31/2012
// ----------------------------------------------------------------------------
template <typename ScalarT>
void akima<ScalarT>::evalDeriv2 (
   const std::vector<ScalarT> & xa,
   const std::vector<ScalarT> & ya,
   const ScalarT & x,
   ScalarT & ypp) const
{
  size_t size = xa.size();
  size_t index = this->binarySearch (xa, x, 0, size - 1);

  const ScalarT delx = x - xa[index];
  ypp = 2.0 * p2[index] + 6.0 * p3[index] * delx;
  return;
}

//-----------------------------------------------------------------------------
// Function      : akima<ScalarT>::evalInteg
// Purpose       :
// Special Notes : 
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 1/31/2012
// ----------------------------------------------------------------------------
template <typename ScalarT>
void akima<ScalarT>::evalInteg (
   const std::vector<ScalarT> & xa,
   const std::vector<ScalarT> & ya,
   const ScalarT & ai,
   const ScalarT & bi,
   ScalarT & result) const
{
  size_t size = xa.size();
  size_t index_a = this->binarySearch (xa, ai, 0, size - 1);
  size_t index_b = this->binarySearch (xa, bi, 0, size - 1);
  result = 0.0;

  // interior intervals
  for(size_t i=index_a; i<=index_b; i++)
  {
    const ScalarT x_hi = xa[i + 1];
    const ScalarT x_lo = xa[i];
    const ScalarT y_lo = ya[i];
    const ScalarT dx = x_hi - x_lo;
    if(dx != 0.0)
    {
      if (i == index_a || i == index_b)
      {
        ScalarT x1 = (i == index_a) ? ai : x_lo;
        ScalarT x2 = (i == index_b) ? bi : x_hi;
        result += this->integ_eval (y_lo, p1[i], p2[i], p3[i], x_lo, x1, x2);
      }
      else
      {
        result += dx * (y_lo + dx*(0.5*p1[i] + dx*(p2[i]/3.0 + 0.25*p3[i]*dx)));
      }
    }
    else
    {
      result = 0.0;
      return;
    }
  }
  return;
}

//-----------------------------------------------------------------------------
// Class         : wodicka spline class
// Purpose       : 
// Special Notes : References:
//
// A variation on Akima splines, which is used in this class, can be 
// found in chapter 13 of the following book:
//
//   "Numerical Algorithms in C" by Gisela Engeln-Mullges and Frank Uhlig
//   Springer 1996.
//   DOI: 10.1007/978-3-642-61074-5
//
//  Two versions of akima splines are presented in the book - one 
//  with rounded corners, and one without.  These versions were 
//  (apparently - it is written in German) described in a 
//  paper by Reinhard Wodicka in 1991:
//
//    WODICKA, R: Ergänzungen zu Akima"s Steigungsformel, Mit teilungen 
//    aus dem Mathem. Seminar Giessen, Heft 203, Selbstver lag des 
//    Mathematischen Instituts, Giessen 1991.
//
// This class implements the non-rounded variant.   It should be functionally
// nearly identical to the "akima" class which based on the original 1970 paper.
// As such this class may be redundant with that one, but I am keeping it around
// for historical reasons.
//
// Creator       : Eric Keiter, SNL
// Creation Date : 9/13/2017
//-----------------------------------------------------------------------------
template <typename ScalarT>
class wodicka: public interpolator<ScalarT>
{
public:
  wodicka () {};
  ~wodicka () { clear(); };

  void init ( const std::vector<ScalarT> & xa, const std::vector<ScalarT> & ya);

  void clear () { p1.clear(); p2.clear(); p3.clear(); m.clear(); };

  void eval (
     const std::vector<ScalarT> & xa, const std::vector<ScalarT> & ya,
     const ScalarT & x, ScalarT & y) const;

  void evalDeriv (
     const std::vector<ScalarT> & xa, const std::vector<ScalarT> & ya,
     const ScalarT & x, ScalarT & dydx) const;

  void evalDeriv2 (
     const std::vector<ScalarT> & xa, const std::vector<ScalarT> & ya,
     const ScalarT & x, ScalarT & ypp) const;

  void evalInteg (
     const std::vector<ScalarT> & xa, const std::vector<ScalarT> & ya,
     const ScalarT & a, const ScalarT & b, ScalarT & result) const;

  void getCoefs (
     const std::vector<ScalarT> & xa, const std::vector<ScalarT> & ya,
        size_t index, std::vector<ScalarT> & coefs)
    {
      coefs.clear();
      if (index < xa.size())
      {
        coefs.resize(4,0.0);
        coefs[0] = ya[index];
        coefs[1] = p1[index];
        coefs[2] = p2[index];
        coefs[3] = p3[index];
      }
    };

public:
  std::vector<ScalarT>  p1;
  std::vector<ScalarT>  p2;
  std::vector<ScalarT>  p3;
  std::vector<ScalarT>  m;
};

//-----------------------------------------------------------------------------
// Function      : wodicka<ScalarT>::init
// Purpose       :
// Special Notes : 
//
// Note that the 'm' indices are offset by 2 on either side of the array,
// compared with the x indices.  ie, m[2] corresponds to x[0], etc.  This is 
// because the Akima algorithm requires 5 points for a given data 
// point - two below and two above.  At each boundary, two additional 
// artificial points which approximate the edge-based derivatives had to be 
// added.
//
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 9/13/2017
// ----------------------------------------------------------------------------
template <typename ScalarT>
void wodicka<ScalarT>::init (
   const std::vector<ScalarT> & xa,
   const std::vector<ScalarT> & ya)
{
  size_t size = xa.size();
  size_t i;

  if (size<=0)
  {
    Report::DevelFatal().in("wodicka<ScalarT>::init") << "Array size  = " << size << ".  Inteprolation failed";
  }

  if (p1.size() != size) p1.resize(size,0.0);
  if (p2.size() != size) p2.resize(size,0.0);
  if (p3.size() != size) p3.resize(size,0.0);

  // m represents edges between points.  So, original size of m is n-1, if the 
  // size of xa and ya is n.  Therefore, when adding the additional 4 points to 
  // m, required by the method, the size of m should be n+3.
  if (m.size() != size+3) 
  {
    m.resize(size+3,0.0);
  }

  // set up the edge derivatives
  if (xa.size() > 1)
  {
    for (i = 0; i < size-1; i++)
    {
      m[i+2] = (ya[i + 1] - ya[i]) / (xa[i + 1] - xa[i]);
    }
  }

  // As described in the original paper, boundary conditions are handled slightly differently.  
  // The slopes for two extra segments on either side of the array need to be estimated.
  //
  // For a non-periodic curve:
  //
  //    indices from 1970 paper:                  indices used here:
  //   m_{-1}  = 3*m_{1} − 2*m_{2}      -> m_{0}   = 3*m_{2} − 2*m_{3}   (adding one)
  //   m_{0}   = 2*m_{1} − m_{2}        -> m_{1}   = 2*m_{1} − m_{2}     (adding one)
  //   m_{n}   = 2*m_{n−1} − m_{n−2}    -> m_{n+1} = 2*m_{n} − m_{n−1}   (adding one)
  //   m_{n+1} = 3*m_{n−1} − 2*m_{n−2}  -> m_{n+2} = 3*m_{n] - 2*m_{n-1} (adding one)
  //
  // For a periodic curve: 
  //
  //    indices from 1970 paper:          indices used here:
  //   m_{−1}  = m_{n−2}        ->  m_0 = m_{n-1}
  //   m_{0}   = m_{n−1}        ->  m_1 = m_{n}
  //   m_{n}   = m_{1}          ->  m_{n+1} = m_{2}
  //   m_{n+1} = m_{2}          ->  m_{n+2} = m_{3}    
  //
  // Assuming non-periodic boundary conditions:
  m[0] = 3.0 * m[2] - 2.0 * m[3];
  m[1] = 2.0 * m[2] - m[3];
  m[size + 1] = 2.0 * m[size] - m[size-1];
  m[size + 2] = 3.0 * m[size] - 2.0 * m[size-1];

  // Periodic boundary conditions:
  //m[0] = m[size-1]
  //m[1] = m[size]
  //m[size + 1] =  m[2]
  //m[size + 2] =  m[3]
 
  ScalarT t1 = 0.0;
  ScalarT t2 = 0.0;

  for (i = 0; i < size-1; i++)
  {
    // Wodicka variation, as documented in the book, "Numerical Algorithms in C" 

    // this is the denominator of eq.(1) in original Akima 1970 paper
    const ScalarT denom = std::abs(m[i + 3] - m[i+2]) + std::abs(m[i + 1] - m[i]);
    if (denom == 0.0)
    {
      // if the denominator of eq.(1) is zero, then this is a special case
      p1[i] = m[i+2];
      p2[i] = 0.0;
      p3[i] = 0.0;
    }
    else
    {
      const ScalarT denom_next = std::abs(m[i + 4] - m[i + 3]) + std::abs(m[i+2] - m[i + 1]);
      const ScalarT alpha_i = std::abs(m[i + 1] - m[i]) / denom;
        
      t1 = (1.0 - alpha_i) * m[i + 1] + alpha_i * m[i+2];

      if (denom_next == 0.0)
      {
        t2 = m[i+2];
      }
      else
      {
        ScalarT alpha_ip1 = std::abs(m[i+2] - m[i + 1]) / denom_next;
        t2 = (1.0 - alpha_ip1) * m[i+2] + alpha_ip1 * m[i + 3];
      }
      
      // see equations 3-6 of the original Akima 1970 paper.  
      // (p0 is used implicitly in the eval function, so not computed/stored)
      const ScalarT dx = xa[i + 1] - xa[i];
      p1[i] = t1;
      p2[i] = (3.0 * m[i+2] - 2.0 * t1 - t2) / dx;
      p3[i] = (t1 + t2 - 2.0 * m[i+2]) / (dx * dx);
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : wodicka<ScalarT>::eval
// Purpose       :
// Special Notes : 
//
//  y(x) = p0+ pl(x-xl) + p2(x-x1)^2+ p3(x-x1)^3
//
//  where:
//
//  p0 = y1 = y(x1)  (result of the binary search)
//  p1 = t1 - estimated slope dydx at x1, based on equation 1 from original paper
//  p2 = [3(y2 - y,)/(x2 - xl) - 2tl - t2]/(x2 - xl)
//  p3 = [t1 + t2 - 2(y2 - yl)/(x2 - xl)]/(x2 - xl)^2
//
//  p0,p1,p2 and p3 were all computed in the "init" function.
//
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 9/13/2017
// ----------------------------------------------------------------------------
template <typename ScalarT>
void wodicka<ScalarT>::eval (
   const std::vector<ScalarT> & xa,
   const std::vector<ScalarT> & ya,
   const ScalarT & x,
   ScalarT & y) const
{
  size_t size = xa.size();
  size_t index = this->binarySearch (xa, x, 0, size - 1);

  const ScalarT delx = x - xa[index];
  y = ya[index] + delx * (p1[index] + delx * (p2[index] + p3[index] * delx));
  return;
}

//-----------------------------------------------------------------------------
// Function      : wodicka<ScalarT>::evalDeriv
// Purpose       :
// Special Notes : 
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 9/13/2017
// ----------------------------------------------------------------------------
template <typename ScalarT>
void wodicka<ScalarT>::evalDeriv (
   const std::vector<ScalarT> & xa,
   const std::vector<ScalarT> & ya,
   const ScalarT & x,
   ScalarT & dydx) const
{
  size_t size = xa.size();
  size_t index = this->binarySearch (xa, x, 0, size - 1);

  ScalarT delx = x - xa[index];
  dydx = p1[index] + delx * (2.0 * p2[index] + 3.0 * p3[index] * delx);
  return;
}

//-----------------------------------------------------------------------------
// Function      : wodicka<ScalarT>::evalDeriv2
// Purpose       :
// Special Notes : 
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 9/13/2017
// ----------------------------------------------------------------------------
template <typename ScalarT>
void wodicka<ScalarT>::evalDeriv2 (
   const std::vector<ScalarT> & xa,
   const std::vector<ScalarT> & ya,
   const ScalarT & x,
   ScalarT & ypp) const
{
  size_t size = xa.size();
  size_t index = this->binarySearch (xa, x, 0, size - 1);

  const ScalarT delx = x - xa[index];
  ypp = 2.0 * p2[index] + 6.0 * p3[index] * delx;
  return;
}

//-----------------------------------------------------------------------------
// Function      : wodicka<ScalarT>::evalInteg
// Purpose       :
// Special Notes : 
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 9/13/2017
// ----------------------------------------------------------------------------
template <typename ScalarT>
void wodicka<ScalarT>::evalInteg (
   const std::vector<ScalarT> & xa,
   const std::vector<ScalarT> & ya,
   const ScalarT & ai,
   const ScalarT & bi,
   ScalarT & result) const
{
  size_t size = xa.size();
  size_t index_a = this->binarySearch (xa, ai, 0, size - 1);
  size_t index_b = this->binarySearch (xa, bi, 0, size - 1);
  result = 0.0;

  // interior intervals
  for(size_t i=index_a; i<=index_b; i++)
  {
    const ScalarT x_hi = xa[i + 1];
    const ScalarT x_lo = xa[i];
    const ScalarT y_lo = ya[i];
    const ScalarT dx = x_hi - x_lo;
    if(dx != 0.0)
    {
      if (i == index_a || i == index_b)
      {
        ScalarT x1 = (i == index_a) ? ai : x_lo;
        ScalarT x2 = (i == index_b) ? bi : x_hi;
        result += this->integ_eval (y_lo, p1[i], p2[i], p3[i], x_lo, x1, x2);
      }
      else
      {
        result += dx * (y_lo + dx*(0.5*p1[i] + dx*(p2[i]/3.0 + 0.25*p3[i]*dx)));
      }
    }
    else
    {
      result = 0.0;
      return;
    }
  }
  return;
}

//-----------------------------------------------------------------------------
// Class         : cubic spline class
// Purpose       :
// Special Notes : 
// Creator       : Eric Keiter, SNL
// Creation Date : 1/31/2012
//-----------------------------------------------------------------------------
template <typename ScalarT>
class cubicSpline: public interpolator<ScalarT>
{
public:
  cubicSpline () {};
  ~cubicSpline () { clear(); };

  void init ( const std::vector<ScalarT> & xa, const std::vector<ScalarT> & ya);

  void clear () { y2.clear(); };

  void eval (
     const std::vector<ScalarT> & xa, const std::vector<ScalarT> & ya,
     const ScalarT & x, ScalarT & y) const;

  void evalDeriv (
     const std::vector<ScalarT> & xa, const std::vector<ScalarT> & ya,
     const ScalarT & x, ScalarT & dydx) const;

  void evalDeriv2 (
     const std::vector<ScalarT> & xa, const std::vector<ScalarT> & ya,
     const ScalarT & x, ScalarT & ypp) const;

  // evalInteg is not implemented
  void getCoefs (
     const std::vector<ScalarT> & xa, const std::vector<ScalarT> & ya,
        size_t index, std::vector<ScalarT> & coefs);
public:
  std::vector<ScalarT> y2;
};

//-----------------------------------------------------------------------------
// Function      : cubicSpline<ScalarT>::init
// Purpose       :
// Special Notes : 
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 1/31/2012
//-----------------------------------------------------------------------------
template <typename ScalarT>
void cubicSpline<ScalarT>::init
(const std::vector<ScalarT> & xa,
 const std::vector<ScalarT> & ya)
{
  if (xa.size()<=0)
  {
    Report::DevelFatal().in("cubicSpline<ScalarT>::init") << "Array size  = " << xa.size() << ".  Inteprolation failed";
  }

  if (y2.size() != xa.size())
  {
    y2.resize(xa.size());
  }

  ScalarT p=0; 
  ScalarT qn=0; 
  ScalarT sig=0; 
  ScalarT un=0;
  int n = y2.size(); 
  std::vector <ScalarT> u(n-1,0.0);

  // natural boundary condition
  y2[0] = 0.0;
  y2[n-1] = 0.0;

  // Tridiagonal solve.
  for (int i=1; i<n-1; i++)
  {
    sig = (xa[i]-xa[i-1])/(xa[i+1]-xa[i-1]);
    p = sig*y2[i-1] + 2.0;
    y2[i] = (sig-1.0)/p;
    u[i] = (ya[i+1]-ya[i])/(xa[i+1]-xa[i]) - (ya[i]-ya[i-1])/(xa[i]-xa[i-1]); // 2nd deriv
    u[i] = (6.0*u[i]/(xa[i+1]-xa[i-1]) - sig*u[i-1])/p;
  }

  for (int l=n-2; l>=0; l--)
  {
    y2[l] = y2[l]*y2[l+1]+u[l];
  }
};

//-----------------------------------------------------------------------------
// Function      : cubicSpline<ScalarT>::eval
// Purpose       :
// Special Notes : 
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 1/31/2012
//-----------------------------------------------------------------------------
template <typename ScalarT>
void cubicSpline<ScalarT>::eval(
   const std::vector<ScalarT> & xa,
   const std::vector<ScalarT> & ya,
   const ScalarT & x_position,
   ScalarT & y_spline) const
{
  int n = xa.size();
  
  ScalarT h = 0.0; 
  ScalarT a = 0.0; 
  ScalarT b = 0.0;

  int k = 0; int klo = 0; int khi = n-1;
  while (khi-klo > 1)
  {
    k = (khi+klo) >> 1;
    if(greaterThan(xa[k],x_position)) khi=k;
    else klo=k;
  }  
  h = xa[khi] - xa[klo];

  if (h == 0.0)
  {
    // if out of range, then use the formula for dy/dx to extrapolate
    // beyond the range.  
    if (khi == 0)
    {
      ScalarT h0 = xa[1]-xa[0];
      ScalarT dx = x_position - xa[0];
      ScalarT dydx = (ya[1]-ya[0])/h0 - h0*y2[0]/3.0 - h0*y2[1]/6.0;
      y_spline = ya[0] + dx * dydx;
    }
    else if (klo == n-1)
    {
      ScalarT h1 = xa[n-1]-xa[n-2];
      ScalarT dx = x_position - xa[n-1];
      ScalarT dydx = (ya[n-1]-ya[n-2])/h1 + h1*y2[n-2]/6.0 + h1*y2[n-1]/3.0;
      y_spline = ya[n-1] + dx * dydx;
    }
  }
  else
  {
    a = (xa[khi] - x_position)/h;
    b = (x_position - xa[klo])/h;
    // cubic spline polynomial: 
    y_spline = a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2[klo] + (b*b*b-b)*y2[khi])*(h*h)/6.0;

  }
}

//-----------------------------------------------------------------------------
// Function      : cubicSpline<ScalarT>::evalDeriv
// Purpose       :
// Special Notes : 
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 1/31/2012
//-----------------------------------------------------------------------------
template <typename ScalarT>
void cubicSpline<ScalarT>::evalDeriv(
   const std::vector<ScalarT> & xa,
   const std::vector<ScalarT> & ya,
   const ScalarT & x_position,
   ScalarT & dydx_spline) const
{
  int n = xa.size();
  // Find the right place in the table by means of bisection.
  ScalarT h = 0.0; ScalarT a = 0.0; ScalarT b = 0.0;
  int k = 0; int klo = 0; int khi = n-1;
  while (khi-klo > 1)
  {
    k = (khi+klo) >> 1;
    if (greaterThan(xa[k], x_position)) khi=k;
    else klo=k;
  }
  h = xa[khi] - xa[klo];
  if (h == 0.0)
  {
    // if out of range, then use the formula for dy/dx to extrapolate
    // beyond the range.  
    if (khi == 0)
    {
      ScalarT h0 = xa[1]-xa[0];
      dydx_spline = (ya[1]-ya[0])/h0 - h0*y2[0]/3.0 - h0*y2[1]/6.0;
    }
    else if (klo == n-1)
    {
      ScalarT h1 = xa[n-1]-xa[n-2];
      dydx_spline = (ya[n-1]-ya[n-2])/h1 + h1*y2[n-2]/6.0 + h1*y2[n-1]/3.0;
    }
  }
  else
  {
    a = (xa[khi] - x_position)/h;
    b = (x_position - xa[klo])/h;

    // derivative:  (formula 3.3.5 from numerical recipies in C)
    dydx_spline =  (ya[khi]-ya[klo])/h - ((3.0*a*a-1.0)*y2[klo] - (3.0*b*b-1.0)*y2[khi])*h/6.0;
  }
}

//-----------------------------------------------------------------------------
// Function      : cubicSpline<ScalarT>::evalDeriv2
// Purpose       :
// Special Notes : 
//
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 1/31/2012
//-----------------------------------------------------------------------------
template <typename ScalarT>
void cubicSpline<ScalarT>::evalDeriv2(
   const std::vector<ScalarT> & xa,
   const std::vector<ScalarT> & ya,
   const ScalarT & x_position,
   ScalarT & ypp) const
{
  int n = xa.size();
  // Find the right place in the table by means of bisection.
  ScalarT h = 0.0; ScalarT a = 0.0; ScalarT b = 0.0;
  int k = 0; int klo = 0; int khi = n-1;
  while (khi-klo > 1)
  {
    k = (khi+klo) >> 1;
    if (greaterThan(xa[k], x_position))
    {
      khi=k;
    }
    else 
    {
      klo=k;
    }
  }
  h = xa[khi] - xa[klo];
  if (h == 0.0)
  {
    // if out of range, assume no curvature.
    if (khi == 0)
    {
      ypp = 0.0;
    }
    else if (klo == n-1)
    {
      ypp = 0.0;
    }
  }
  else
  {
    a = (xa[khi] - x_position)/h;
    b = (x_position - xa[klo])/h;

    ypp =  a*y2[klo] + b*y2[khi];
  }
}

//-----------------------------------------------------------------------------
// Function      : cubicSpline<ScalarT>::getCoefs
// Purpose       : provides the cubic polynomial coefficients for a given cell.
//
// Special Notes : the cubic spline evaluation formula is not in the standard
//                 polynomial form. So, do re-arrangements to get it there.
//
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 10/03/2024
//-----------------------------------------------------------------------------
template <typename ScalarT>
void cubicSpline<ScalarT>::getCoefs(
     const std::vector<ScalarT> & xa, const std::vector<ScalarT> & ya,
    size_t index, std::vector<ScalarT> & coefs)
{
  int n = xa.size();

  if (index < n)
  {
    coefs.resize(4,0.0);
  }

  ScalarT h = 0.0;
  ScalarT a = 0.0;
  ScalarT b = 0.0;

  int k = 0;
  int klo = index; 
  int khi = index+1;

  if (klo < 0) { klo = 0; }
  if (khi < 0) { khi = 0; }

  if (klo > n-1) { klo = n-1; }
  if (khi > n-1) { khi = n-1; }

  h = xa[khi] - xa[klo];

  if (h == 0.0)
  {
    // if out of range, then use the formula for dy/dx to extrapolate
    // beyond the range.
    if (khi == 0)
    {
#if 0
      ScalarT h0 = xa[1]-xa[0];
      ScalarT dx = x_position - xa[0];
      y_spline = ya[0] + dx * dydx;
#endif
      // y_spline = ya[0] + dx * dydx = ya[0] + (x-xa[0]) * dydx

      ScalarT h0 = xa[1]-xa[0];
      ScalarT dydx = (ya[1]-ya[0])/h0 - h0*y2[0]/3.0 - h0*y2[1]/6.0;

      ScalarT p0 = ya[0] - xa[0] * dydx;
      ScalarT p1 = dydx;
      ScalarT p2 = 0.0;
      ScalarT p3 = 0.0;
    }
    else if (klo == n-1)
    {
#if 0
      ScalarT h1 = xa[n-1]-xa[n-2];
      ScalarT dx = x_position - xa[n-1];
      ScalarT dydx = (ya[n-1]-ya[n-2])/h1 + h1*y2[n-2]/6.0 + h1*y2[n-1]/3.0;
      y_spline = ya[n-1] + dx * dydx;
#endif

      ScalarT h1 = xa[n-1]-xa[n-2];
      ScalarT dydx = (ya[n-1]-ya[n-2])/h1 + h1*y2[n-2]/6.0 + h1*y2[n-1]/3.0;

      ScalarT p0 = ya[n-1] - xa[n-1]*dydx;
      ScalarT p1 = dydx;
      ScalarT p2 = 0.0;
      ScalarT p3 = 0.0;

      coefs[0] = p0;
      coefs[1] = p1;
      coefs[2] = p2;
      coefs[3] = p3;
    }
  }
  else
  {
#if 0
    a = (xa[khi] - x_position)/h;
    b = (x_position - xa[klo])/h;
    // cubic spline polynomial: 
    y_spline = a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2[klo] + (b*b*b-b)*y2[khi])*(h*h)/6.0;
#endif

    if (index<n)
    {
      // the cubic spline evaluation formula is not in the standard polynomial form.
      // So, do re-arrangements to get it there.
      //a = (xa[khi] - x_position)/h;
      //b = (x_position - xa[klo])/h;

      ScalarT p0 = 0.0;
      ScalarT p1 = 0.0;
      ScalarT p2 = 0.0;
      ScalarT p3 = 0.0;

      p0 += (xa[khi] )/h * ya[klo];   // from the a*ya[klo] term (1st term)
      p0 += (- xa[klo])/h * ya[khi];  // from the b*ya[khi] term (2nd term)
      p0 += + y2[klo]/h * std::pow(xa[khi], 3.0) / 6.0 - y2[klo] * h * xa[khi] / 6.0;  // from the 3rd term: ((a*a*a-a)*y2[klo])*(h*h)/6.0
      p0 += - y2[khi]/h * std::pow(xa[klo], 3.0) / 6.0 + y2[khi] * h * xa[klo] / 6.0;  // from the 4th term: ((b*b*b-b)*y2[khi])*(h*h)/6.0

      p1 += ( - 1.0)/h * ya[klo];
      p1 += (1.0)/h * ya[khi];
      p1 += + (-y2[klo] / h * xa[khi] * xa[khi] / 2.0 + y2[klo] * h / 6.0);
      p1 += + (y2[khi] / h * xa[klo] * xa[klo] / 2.0 - y2[khi] * h / 6.0);

      p2 += + y2[klo] / h * xa[khi] / 2.0;
      p2 += - y2[khi] / h * xa[klo] / 2.0;

      p3 += -y2[klo] / h / 6.0;
      p3 += y2[khi] / h / 6.0;

      coefs[0] = p0;
      coefs[1] = p1;
      coefs[2] = p2;
      coefs[3] = p3;
    }
  }
}


//-----------------------------------------------------------------------------
// Class         : linear interpolation class
// Purpose       :
// Special Notes : 
// Creator       : Eric Keiter, SNL
// Creation Date : 1/31/2012
//-----------------------------------------------------------------------------
template <typename ScalarT>
class linear: public interpolator<ScalarT>
{
public:
  linear () {};

  void init ( const std::vector<ScalarT> & xa, const std::vector<ScalarT> & ya) {};

  void clear () { };

  void eval (
     const std::vector<ScalarT> & xa, const std::vector<ScalarT> & ya,
     const ScalarT & x, ScalarT & y) const;

  void evalDeriv (
     const std::vector<ScalarT> & xa, const std::vector<ScalarT> & ya,
     const ScalarT & x, ScalarT & dydx) const;

  void evalDeriv2 (
     const std::vector<ScalarT> & xa, const std::vector<ScalarT> & ya,
     const ScalarT & x, ScalarT & ypp) const;

  void evalInteg (
     const std::vector<ScalarT> & xa, const std::vector<ScalarT> & ya,
     const ScalarT & a, const ScalarT & b, ScalarT & result) const;

public:

};


//-----------------------------------------------------------------------------
// Function      : linear<ScalarT>::eval
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 1/31/2012
//-----------------------------------------------------------------------------
template <typename ScalarT>
void linear<ScalarT>::eval(
   const std::vector<ScalarT> & xa,
   const std::vector<ScalarT> & ya,
   const ScalarT & x,
   ScalarT & y) const
{
  int n = xa.size();

  // Find the right place in the table by means of bisection.
  ScalarT h = 0.0; ScalarT a = 0.0; ScalarT b = 0.0;
  int k = 0; int klo = 0; int khi = n-1;
  while (khi-klo > 1)
  {
    k = (khi+klo) >> 1;

    if(greaterThan(xa[k],x))
    {
      khi=k;
    }
    else 
    {
      klo=k;
    }
  }
  h = xa[khi] - xa[klo];

  if (h == 0.0)
  {
    if (khi == 0)
    {
      y = xa[khi];
    }
    else if (klo == n-1)
    {
      y = xa[klo];
    }
  }
  else
  {
    ScalarT dx = x - xa[klo];
    ScalarT ya0 = ya[khi] - ya[klo];
    y = (dx/h) * ya0 + ya[klo];
  }
}

//-----------------------------------------------------------------------------
// Function      : linear<ScalarT>::evalDeriv
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 1/31/2012
//-----------------------------------------------------------------------------
template <typename ScalarT>
void linear<ScalarT>::evalDeriv(
   const std::vector<ScalarT> & xa,
   const std::vector<ScalarT> & ya,
   const ScalarT & x,
   ScalarT & dydx) const
{
  int n = xa.size();

  // Find the right place in the table by means of bisection.
  ScalarT h = 0.0; ScalarT a = 0.0; ScalarT b = 0.0;
  int k = 0; int klo = 0; int khi = n-1;
  while (khi-klo > 1)
  {
    k = (khi+klo) >> 1;

    if(greaterThan(xa[k],x))
    {
      khi=k;
    }
    else 
    {
      klo=k;
    }
  }
  h = xa[khi] - xa[klo];

  if (h == 0.0)
  {
    if (khi == 0)
    {
      dydx = 0.0;
    }
    else if (klo == n-1)
    {
      dydx = 0.0;
    }
  }
  else
  {
    ScalarT dx = xa[khi] - xa[klo];
    ScalarT dy  = ya[khi] - ya[klo];
    dydx = dy/dx;
  }
}

//-----------------------------------------------------------------------------
// Function      : linear<ScalarT>::evalDeriv2
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 1/31/2012
//-----------------------------------------------------------------------------
template <typename ScalarT>
void linear<ScalarT>::evalDeriv2(
   const std::vector<ScalarT> & xa,
   const std::vector<ScalarT> & ya,
   const ScalarT & x,
   ScalarT & ypp) const
{
  ypp = 0.0;
}

//-----------------------------------------------------------------------------
// Function      : linear<ScalarT>::evalInteg
// Purpose       :
// Special Notes : 
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 1/31/2012
// ----------------------------------------------------------------------------
template <typename ScalarT>
void linear<ScalarT>::evalInteg (
   const std::vector<ScalarT> & xa,
   const std::vector<ScalarT> & ya,
   const ScalarT & a,
   const ScalarT & b,
   ScalarT & result) const
{

  int size = xa.size();
  int index_a = this->binarySearch (xa, a, 0, size - 1);
  int index_b = this->binarySearch (xa, b, 0, size - 1);

  // endpoints span more than one interval
  result = 0.0;

  // interior intervals
  for(int i=index_a; i<=index_b; i++)
  {
    const ScalarT x_hi = xa[i + 1];
    const ScalarT x_lo = xa[i];
    const ScalarT y_lo = ya[i];
    const ScalarT y_hi = ya[i + 1];
    const ScalarT dx = x_hi - x_lo;

    if(dx != 0.0)
    {
      if (i == index_a || i == index_b)
      {
        ScalarT x1 = (i == index_a) ? a : x_lo;
        ScalarT x2 = (i == index_b) ? b : x_hi;
        const ScalarT D = (y_hi-y_lo)/dx;
        result += (x2-x1) * (y_lo + 0.5*D*((x2-x_lo)+(x1-x_lo)));
      }
      else
      {
        result += 0.5 * dx * (y_lo + y_hi);
      }
    }
  }
  return;
}

//-----------------------------------------------------------------------------
// Class         : quadratic spline interpolation class
// Purpose       : 
// Special Notes : 
// Creator       : Eric Keiter, SNL
// Creation Date : 9/9/2017
//-----------------------------------------------------------------------------
template <typename ScalarT>
class quadSpline: public interpolator<ScalarT>
{
public:
  quadSpline () {};

  void init ( const std::vector<ScalarT> & xa, const std::vector<ScalarT> & ya);

  void clear () { };

  void eval (
     const std::vector<ScalarT> & xa, const std::vector<ScalarT> & ya,
     const ScalarT & x, ScalarT & y) const;

  void evalDeriv (
     const std::vector<ScalarT> & xa, const std::vector<ScalarT> & ya,
     const ScalarT & x, ScalarT & dydx) const;

  void evalDeriv2 (
     const std::vector<ScalarT> & xa, const std::vector<ScalarT> & ya,
     const ScalarT & x, ScalarT & ypp) const;

  // evalInteg not implemented

public:
  std::vector<ScalarT> b;
  std::vector<ScalarT> c;
};

//-----------------------------------------------------------------------------
// Function      : quadSpline<ScalarT>::init
// Purpose       :
// Special Notes : 
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 9/9/2017
//-----------------------------------------------------------------------------
template <typename ScalarT>
void quadSpline<ScalarT>::init
(const std::vector<ScalarT> & xa,
 const std::vector<ScalarT> & ya)
{
  int size = xa.size();

  if (size<=0)
  {
    Report::DevelFatal().in("quadSpline<ScalarT>::init") << "Array size  = " << size << ".  Inteprolation failed";
  }

  if (b.size() != size) b.resize(size,0.0);
  if (c.size() != size) c.resize(size,0.0);

  std::vector<ScalarT> h(size,0.0);
  std::vector<ScalarT> p(size,0.0);

  for(int i=0; i<size-1; i++)
  {
    h[i]=xa[i+1]-xa[i];
    p[i]=(ya[i+1]-ya[i])/h[i];
  }
  c[0]=0;

  //forward:
  for(int i=0; i<size-2; i++)
  {
    c[i+1]=(p[i+1]-p[i]-c[i]*h[i])/h[i+1];
  }

  c[size-2]/=2;

  //bacwkard
  for(int i=size-3; i>=0; i--)
  {
    c[i]=(p[i+1]-p[i]-c[i+1]*h[i+1])/h[i];
  }

  for(int i=0; i<size-1; i++)
  {
    b[i]=p[i]-c[i]*h[i];
  }

  return;
}

//-----------------------------------------------------------------------------
// Function      : quadSpline<ScalarT>::eval
// Purpose       :
// Special Notes : 
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 9/9/2017
//-----------------------------------------------------------------------------
template <typename ScalarT>
void quadSpline<ScalarT>::eval(
   const std::vector<ScalarT> & xa,
   const std::vector<ScalarT> & ya,
   const ScalarT & x_position,
   ScalarT & y_spline) const
{
  int size = xa.size();
  size_t index = this->binarySearch (xa, x_position, 0, size - 1);
  ScalarT h=x_position-xa[index];

  y_spline = ya[index]+h*(b[index]+h*c[index]);
}

//-----------------------------------------------------------------------------
// Function      : quadSpline<ScalarT>::evalDeriv
// Purpose       :
// Special Notes : 
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 9/9/2017
//-----------------------------------------------------------------------------
template <typename ScalarT>
void quadSpline<ScalarT>::evalDeriv(
   const std::vector<ScalarT> & xa,
   const std::vector<ScalarT> & ya,
   const ScalarT & x_position,
   ScalarT & dydx) const
{
  int size = xa.size();
  size_t index = this->binarySearch (xa, x_position, 0, size - 1);
  ScalarT h=x_position-xa[index];

  dydx = b[index]+2*h*c[index];
}

//-----------------------------------------------------------------------------
// Function      : quadSpline<ScalarT>::evalDeriv2
// Purpose       :
// Special Notes : 
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 9/9/2017
//-----------------------------------------------------------------------------
template <typename ScalarT>
void quadSpline<ScalarT>::evalDeriv2(
   const std::vector<ScalarT> & xa,
   const std::vector<ScalarT> & ya,
   const ScalarT & x_position,
   ScalarT & ypp) const
{
  int size = xa.size();
  size_t index = this->binarySearch (xa, x_position, 0, size - 1);
  ypp = 2*c[index];
}

//-----------------------------------------------------------------------------
// Class         : barycentricLagrange interpolation class
// Purpose       : 
// Special Notes : Implements Lagrange interpolation, which has been modified 
// to be "Barycentric".  Based on the method described in:
//
//   "Barycentric Lagrange Interpolation" by Jean-Paul Berrut and 
//   Lloyd N. Trefethen. SIAM REVIEW.  Vol. 46, No. 3, pp. 501–517, 2004.
//
// Note that this method works best if the points used result in a 
// well-conditioned polynomial.  Equally spaced points will not give this result.
// Instead, it is best to choose points schemes which cluster points near the 
// end points of the interval.  As noted in the above paper:
//
// "As is well known in approximation theory, the right approach is to use 
// point sets that are clus- tered at the endpoints of the interval with 
// an asymptotic density proportional to (1 − x2)−1/2 as n → ∞."
//
// And, later in the paper:
//
// "The simplest examples of clustered point sets are the families 
// of Chebyshev points, obtained by projecting equally spaced points 
// on the unit circle down to the unit interval [−1, 1]."
//
// For example, if using Chebyshev points of the second kind, all the 
// coefficients wind up being 1 or -1, except at the end points, where they 
// are 1/2.
//
// Creator       : Eric Keiter, SNL
// Creation Date : 9/13/2017
//-----------------------------------------------------------------------------
template <typename ScalarT>
class barycentricLagrange: public interpolator<ScalarT>
{
public:
  barycentricLagrange () {};
  ~barycentricLagrange () {clear();};

  void init ( const std::vector<ScalarT> & xa, const std::vector<ScalarT> & ya);

  void clear () { w.clear(); }

  void eval (
     const std::vector<ScalarT> & xa, const std::vector<ScalarT> & ya,
     const ScalarT & x, ScalarT & y) const;

  void evalDeriv (
     const std::vector<ScalarT> & xa, const std::vector<ScalarT> & ya,
     const ScalarT & x, ScalarT & dydx) const;

  // evalDeriv2 and evalInteg not implemented

public:
  std::vector<ScalarT>  w;
};

//-----------------------------------------------------------------------------
// Function      : barycentricLagrange<ScalarT>::init
// Purpose       :
// Special Notes : 
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 9/13/2017
// ----------------------------------------------------------------------------
template <typename ScalarT>
void barycentricLagrange<ScalarT>::init (
   const std::vector<ScalarT> & xa,
   const std::vector<ScalarT> & ya)
{
  size_t size = xa.size();

  if (size<=0)
  {
    Report::DevelFatal().in("barycentricLagrange<ScalarT>::init") << "Array size  = " << size << ".  Inteprolation failed";
  }

  w.resize(size,0.0);

  // Compute the weights using formula 3.2 from the paper
  for (size_t j=0;j<size;++j)
  {
    w[j] = 1.0;
    for (size_t k=0;k<size;++k)
    {
      if (!(j==k))
      {
        w[j] *= (xa[j]-xa[k]);
      }
    }

    w[j] = 1.0/w[j];
  }
}

//-----------------------------------------------------------------------------
// Function      : barycentricLagrange<ScalarT>::eval
// Purpose       :
// Special Notes : 
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 9/13/2017
// ----------------------------------------------------------------------------
template <typename ScalarT>
void barycentricLagrange<ScalarT>::eval (
   const std::vector<ScalarT> & xa,
   const std::vector<ScalarT> & ya,
   const ScalarT & x,
   ScalarT & y) const
{
  size_t size = xa.size();

  // Compute l(x) from equation 3.1 of the paper:
  ScalarT l = 1.0;
  int foundExact = -1;
  for (size_t j=0;j<size;++j)
  {
    ScalarT dx = x-xa[j];
    l *= dx;
    if (dx == 0.0)
    {
      foundExact = j;
      break;
    }
  }

  // if dx == 0, then algorithm should return precisely one of the
  // original dataset points.  When dx is zero, the formula doesn't work
  // as it returns nan.
  if (foundExact >= 0)
  {
    y = ya[foundExact];
  }
  else
  {
    // Compute p(x) from equation 3.3 of the paper.  (Here this is y, not p)
    y = 0.0;
    for (size_t j=0;j<size;++j)
    {
      y += (w[j]/(x-xa[j])) * ya[j];
    }
    y *= l;
  }

  return;
}

//-----------------------------------------------------------------------------
// Function      : barycentricLagrange<ScalarT>::evalDeriv
// Purpose       :
// Special Notes : 
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 9/13/2017
// ----------------------------------------------------------------------------
template <typename ScalarT>
void barycentricLagrange<ScalarT>::evalDeriv (
   const std::vector<ScalarT> & xa,
   const std::vector<ScalarT> & ya,
   const ScalarT & x,
   ScalarT & dydx) const
{
  size_t size = xa.size();

  // Compute l(x) from equation 3.1 of the paper:
  ScalarT l = 1.0;
  for (size_t j=0;j<size;++j)
  {
    ScalarT dx = x-xa[j];
    l *= dx;
  }

  // Compute d/dx of l(x) from equation 3.1 of the paper:
  ScalarT dldx = 0.0;
  for (size_t j=0;j<size;++j)
  {
    dldx += l/(x-xa[j]);
  }

  // Compute ddx(p(x)) from equation 3.3 of the paper.  (Here this is y, not p)
  //    interpolated y = l * sum, where the sum is defined in equation 3.3, and l was computed above.
  //   So, derivative dydx = dldx * sum + l * dsumdx
  //
  //   l and dldx were computed above, so now need sum and dsumdx
  ScalarT sum = 0.0;
  ScalarT dsumdx = 0.0;
  for (size_t j=0;j<size;++j)
  {
    ScalarT dx = x-xa[j];
    sum += (w[j]/dx) * ya[j];
    dsumdx -= (w[j]/(dx*dx)) * ya[j];
  }
  dydx = dldx * sum + l * dsumdx;

  return;
}

//-----------------------------------------------------------------------------
// Class         : steffen spline class
// Purpose       :
// Special Notes : References:
//
// The original reference for this type of spline is in this paper:
//
// M.Steffen, "A simple method for monotonic interpolation in one dimension",
// Astron. Astrophys. 239, 443-450 (1990).
//
// Creator       : Eric Keiter, SNL
// Creation Date : 9/25/2024
//-----------------------------------------------------------------------------
template <typename ScalarT>
class steffen: public interpolator<ScalarT>
{
public:
  steffen () {};
  ~steffen () { clear(); };

  ScalarT copysign(const ScalarT x, const ScalarT y)
  {
    if ((x < 0 && y > 0) || (x > 0 && y < 0)) { return -x; }
    return x;
  };

  void init ( const std::vector<ScalarT> & xa, const std::vector<ScalarT> & ya);

  void clear () { a.clear(); b.clear(); c.clear(); d.clear(); yp.clear(); };

  void eval (
     const std::vector<ScalarT> & xa, const std::vector<ScalarT> & ya,
     const ScalarT & x, ScalarT & y) const;

  void evalDeriv (
     const std::vector<ScalarT> & xa, const std::vector<ScalarT> & ya,
     const ScalarT & x, ScalarT & dydx) const;

  void evalDeriv2 (
     const std::vector<ScalarT> & xa, const std::vector<ScalarT> & ya,
     const ScalarT & x, ScalarT & ypp) const;

  void evalInteg (
     const std::vector<ScalarT> & xa, const std::vector<ScalarT> & ya,
     const ScalarT & a, const ScalarT & b, ScalarT & result) const;

  void getCoefs (
     const std::vector<ScalarT> & xa, const std::vector<ScalarT> & ya,
        size_t index, std::vector<ScalarT> & coefs)
    {
      coefs.clear();
      if (index < xa.size())
      {
        // Assume this form/order:  y(x) = coef_0 + coef_l(x-xl) + coef_2(x-x1)^2+ coef_3(x-x1)^3 
        // Given that:  y = a*delx*delx*delx + b*delx*delx + c*delx + d, the order here is d,c,b,a
        coefs.resize(4,0.0);
        coefs[0] = d[index];
        coefs[1] = c[index];
        coefs[2] = b[index];
        coefs[3] = a[index];
      }
    };

public:
  std::vector<ScalarT>  a;
  std::vector<ScalarT>  b;
  std::vector<ScalarT>  c;
  std::vector<ScalarT>  d;
  std::vector<ScalarT>  yp;
};

template <typename ScalarT>
inline ScalarT Xycemax ( ScalarT f1, ScalarT f2) { return f1 > f2 ? f1 : f2; }

template <typename ScalarT>
inline ScalarT Xycemin ( ScalarT f1, ScalarT f2) { return f1 < f2 ? f1 : f2; }

template <typename ScalarT>
inline ScalarT Xyceabs ( ScalarT f1 )  { return (f1 < 0.0)?(-1.0*f1):(f1); }

//-----------------------------------------------------------------------------
// Function      : steffen<ScalarT>::init
// Purpose       :
// Special Notes :
//
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 1/31/2024
// ----------------------------------------------------------------------------
template <typename ScalarT>
void steffen<ScalarT>::init (
   const std::vector<ScalarT> & xa,
   const std::vector<ScalarT> & ya)
{
  size_t size = xa.size();
  size_t i;

  if (a.size() != size) a.resize(size,0.0);
  if (b.size() != size) b.resize(size,0.0);
  if (c.size() != size) c.resize(size,0.0);
  if (d.size() != size) d.resize(size,0.0);
  if (yp.size() != size) yp.resize(size,0.0);

  // Assign interval and slopes for left boundary.
  ScalarT h0 = (xa[1] - xa[0]);
  ScalarT s0 = (ya[1] - ya[0]) / h0;

  yp[0] = s0;

  for (i = 1; i < (size - 1); i++)
  {
      ScalarT pi;

      // equation 6 in paper
      ScalarT hi = (xa[i+1] - xa[i]);
      ScalarT him1 = (xa[i] - xa[i - 1]);

      // equation 7 in paper
      ScalarT si = (ya[i+1] - ya[i]) / hi;
      ScalarT sim1 = (ya[i] - ya[i - 1]) / him1;

      // equation 8 in paper
      pi = (sim1*hi + si*him1) / (him1 + hi);

      // equation 11 in paper
      yp[i] = (copysign(1.0,sim1) + copysign(1.0,si)) * Xycemin(Xyceabs(sim1),Xycemin(Xyceabs(si),Xyceabs(pi)));
  }

  // Assign y' for rightmost boundary
  yp[size-1] = (ya[size - 1] - ya[size - 2]) /
                    (xa[size - 1] - xa[size - 2]);

  for (i = 0; i < (size - 1); i++)
    {
      ScalarT hi = (xa[i+1] - xa[i]);
      ScalarT si = (ya[i+1] - ya[i]) / hi;

      // Equations 2-5 in paper.
      a[i] = (yp[i] + yp[i+1] - 2*si) / hi / hi;
      b[i] = (3*si - 2*yp[i] - yp[i+1]) / hi;
      c[i] = yp[i];
      d[i] = ya[i];
    }

  return;
}


//-----------------------------------------------------------------------------
// Function      : steffen<ScalarT>::eval
// Purpose       :
// Special Notes :
//
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 1/31/2024
// ----------------------------------------------------------------------------
template <typename ScalarT>
void steffen<ScalarT>::eval (
   const std::vector<ScalarT> & xa,
   const std::vector<ScalarT> & ya,
   const ScalarT & x,
   ScalarT & y) const
{
  size_t size = xa.size();
  size_t index = this->binarySearch (xa, x, 0, size - 1);
  const ScalarT x_lo = xa[index];
  const ScalarT delx = x - x_lo;

  // Horner's scheme
  // y = a*delx*delx*delx + b*delx*delx + c*delx + d;
  y = d[index] + delx*(c[index] + delx*(b[index] + delx*a[index]));

  return;
}


//-----------------------------------------------------------------------------
// Function      : steffen<ScalarT>::evalDeriv
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 1/31/2024
// ----------------------------------------------------------------------------
template <typename ScalarT>
void steffen<ScalarT>::evalDeriv (
   const std::vector<ScalarT> & xa,
   const std::vector<ScalarT> & ya,
   const ScalarT & x,
   ScalarT & dydx) const
{
  size_t size = xa.size();
  size_t index = this->binarySearch (xa, x, 0, size - 1);
  ScalarT x_lo = xa[index];
  ScalarT delx = x - x_lo;
  //ScalarT d = state->d[index];
  // dydx = 3*a*delx*delx*delx + 2*b*delx + c;
  dydx = c[index] + delx*(2*b[index] + delx*3*a[index]);
  return;
}

//-----------------------------------------------------------------------------
// Function      : steffen<ScalarT>::evalDeriv2
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 1/31/2024
// ----------------------------------------------------------------------------
template <typename ScalarT>
void steffen<ScalarT>::evalDeriv2 (
   const std::vector<ScalarT> & xa,
   const std::vector<ScalarT> & ya,
   const ScalarT & x,
   ScalarT & ypp) const
{
  size_t size = xa.size();
  size_t index = this->binarySearch (xa, x, 0, size - 1);
  const ScalarT x_lo = xa[index];
  const ScalarT delx = x - x_lo;
  ypp = 6*a[index]*delx + 2*b[index];
  return;
}

//-----------------------------------------------------------------------------
// Function      : steffen<ScalarT>::evalInteg
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 1/31/2024
// ----------------------------------------------------------------------------
template <typename ScalarT>
void steffen<ScalarT>::evalInteg (
   const std::vector<ScalarT> & xa,
   const std::vector<ScalarT> & ya,
   const ScalarT & ai,
   const ScalarT & bi,
   ScalarT & result) const
{
  // ai and bi are the boundaries of the integration.
  size_t i, index_a, index_b;
  size_t size = xa.size();

  // Find the data points in the xa that are nearest to the desired
  // a and b integration boundaries.
  index_a = this->binarySearch (xa, ai, 0, size - 1);
  index_b = this->binarySearch (xa, bi, 0, size - 1);

  result = 0.0;

  for(i=index_a; i<=index_b; i++)
  {
    const ScalarT x_hi = xa[i + 1];
    const ScalarT x_lo = xa[i];
    const ScalarT dx = x_hi - x_lo;
    if(dx != 0.0)
    {
      // check if we are at a boundary point, so take the
      // a and b parameters instead of the data points.
      ScalarT zero = 0.0;
      ScalarT x1 = (i == index_a) ? (ai-x_lo) : zero;
      ScalarT x2 = (i == index_b) ? (bi-x_lo) : (x_hi-x_lo);

      result += (1.0/4.0)*a[i]*(x2*x2*x2*x2 - x1*x1*x1*x1)
                +(1.0/3.0)*b[i]*(x2*x2*x2 - x1*x1*x1)
                +(1.0/2.0)*c[i]*(x2*x2 - x1*x1)
                +d[i]*(x2-x1);
    }
    else
    {
      result = 0.0;
      return;
    }
  }
  return;
}

} // namespace Util
} // namespace Xyce

#endif


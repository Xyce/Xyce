//-------------------------------------------------------------------------
//   Copyright 2002-2020 National Technology & Engineering Solutions of
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


#include <Xyce_config.h>

#include <N_ERH_ErrorMgr.h>
#include <N_ANP_UQSupport.h>
#include <N_ANP_StepEvent.h>
#include <N_LOA_CktLoader.h>
#include <N_UTL_FeatureTest.h>
#include <N_IO_CmdParse.h>
#include <N_NLS_SensitivityResiduals.h>

#include <N_PDS_Serial.h>
#include <N_PDS_MPI.h>

#include <Teuchos_BLAS.hpp>
#include <Teuchos_Utils.hpp>
#include <Teuchos_LAPACK.hpp>

#if __cplusplus>=201103L
// note, this only works with C++11!
#include <random>
#else
#include <N_UTL_RandomNumbers.h>
#endif

#if( defined HAVE__ISNAN_AND__FINITE_SUPPORT )
#include <float.h>
#define isnan(x) _isnan(x)
#define isinf(x) (!_finite(x))
#else
#define isnan(x) std::isnan(x)
#define isinf(x) std::isinf(x)
#endif

namespace Xyce {
namespace Analysis {
namespace UQ {

//-------------------------------------------------------------------------------
// Function      : computeStats
// Purpose       : Computes a bunch of statistical moments from sample data
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 8/1/2018
//-----------------------------------------------------------------------------
void computeStats (const std::vector<double> & values, statisticalMoments & sm)
{
  double sum=0.0;
  double sumSq=0.0;
  double sumCub=0.0;
  double sum4=0.0;

  sm.max = -Xyce::Util::MachineDependentParams::DoubleMax();
  sm.min = +Xyce::Util::MachineDependentParams::DoubleMax();

  int numSamples = values.size();
  if (numSamples==0) return;

  for (int i=0;i<numSamples;++i)
  {
    double val = values[i];
    sum += val;
    sumSq += val*val;
    sumCub += val*val*val;
    sum4 += val*val*val*val;

    if (sm.max < val) sm.max = val;
    if (sm.min > val) sm.min = val;
  }

  double mean = sum / static_cast<double>(numSamples);
  double meanSq = sumSq / static_cast<double>(numSamples);
  double meanCub = sumCub / static_cast<double>(numSamples);
  double mean4 = sum4 / static_cast<double>(numSamples);

  double variance = 0.0;
  //double variance = meanSq - mean*mean; 
  // variance should never be neagive, but this (above) expression occasionally is.
  // Compute it in a slightly more expensive way that ensures it is positive.
  for (int i=0;i<numSamples;++i) { variance += (values[i]-mean)*(values[i]-mean); }
  variance /= static_cast<double>(numSamples);
  double stddev = std::sqrt(variance);

  double mu3=meanCub-3*mean*meanSq+2*mean*mean*mean;
  double sigma3=stddev*stddev*stddev;
  double mu4=mean4-4*mean*meanCub+6*mean*mean*meanSq-3*mean*mean*mean*mean;
  double sigma4=stddev*stddev*stddev*stddev;

  double skew = mu3/sigma3;
  double kurtosis = mu4/sigma4;

  // hopefully these conditions are never triggered.  The main 
  // use case in which this has been observed to happen is if 
  // variance is negative, which produces stddev=nan.  That problem
  // has been mitigated, above.
  if ( isinf(mean) || isnan(mean) ) { mean = 0.0; }
  if ( isinf(stddev) || isnan(stddev) ) { stddev = 0.0; }
  if ( isinf(variance) || isnan(variance) ) { variance = 0.0; }
  if ( isinf(skew) || isnan(skew) ) { skew = 0.0; }
  if ( isinf(kurtosis) || isnan(kurtosis) ) { kurtosis = 0.0; }

  sm.mean = mean;
  sm.stddev = stddev;
  sm.variance = variance;
  sm.skew = skew;
  sm.kurtosis = kurtosis;
}

#if Xyce_STOKHOS_ENABLE
//-------------------------------------------------------------------------------
// Function      : solveRegressionPCE
// Purpose       : solves for the PCE coefficients using regression for a single output 
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 8/1/2018
//-------------------------------------------------------------------------------
void solveRegressionPCE(
    const int numParams, // dimension (number of params)
    const std::vector< std::vector<double> > & x,
    const std::vector<double> & f,
    Stokhos::OrthogPolyApprox<int,double> & samplePCE
    )
{
  const int N = f.size(); // total number of samples

  using Teuchos::Array;
  using Teuchos::RCP;
  using Teuchos::rcp;

  // To create the vandermonde matrix, I need to be able to evaluate each *term* at x.
  // The linear system to be solved via least squares is:
  //
  //  A * c = f
  //
  //  where c is a vector of coefficients, and f is the result of a function evaluation at a point.
  //
  //  The rectangular matrix A is organized  so that the columns are of size N = sampling_size and
  //  the rows are of size K = # coefficients.
  //
  //  So each row corresponds to a given point (x, or p, or whatever).  Each col corresponds to a
  //  coefficient, or specific polynomial term.
  int basisSize = samplePCE.basis()->size();
  Teuchos::SerialDenseMatrix< int, double > Vmatrix(N,basisSize);
  Vmatrix.putScalar(0.0);

  for( std::vector<double>::size_type ii=0; ii<N; ii++) 
  {
    // set up the point
    Array<double> point(numParams);

    for (int id=0;id<numParams;id++)
    {
      point[id] = x[id][ii];
    }

    Teuchos::Array<double> basis_vals(basisSize);
    samplePCE.basis()->evaluateBases(point,basis_vals);

    for( int jj=0; jj<basisSize; jj++) 
    {
      Vmatrix(ii,jj) = basis_vals[jj];
    }
  }

  // Update so RHS is populated with sampling points results.
  Teuchos::SerialDenseMatrix< int, double > RHS(N,1);
  RHS.putScalar(0.0);

  for (int irhs=0;irhs<N;irhs++)
  {
    RHS(irhs,0) = f[irhs];
  }

  // Solve least squares problem using LAPACK
  Teuchos::LAPACK< int, double > lapack;

  double lworkScalar(1.0);
  int info = 0;
  lapack.GELS('N', Vmatrix.numRows(), Vmatrix.numCols(), RHS.numCols(),
                   Vmatrix.values(),  Vmatrix.numRows(), RHS.values(),    RHS.numRows(),
                   &lworkScalar, -1, &info);

  TEUCHOS_TEST_FOR_EXCEPTION(info != 0, std::logic_error,
                             "_GELSS workspace query returned INFO = "
                             << info << " != 0.");

  const int lwork = static_cast<int> (lworkScalar);

  TEUCHOS_TEST_FOR_EXCEPTION(lwork < 0, std::logic_error,
                             "_GELSS workspace query returned LWORK = "
                             << lwork << " < 0.");

  // Allocate workspace.  Size > 0 means &work[0] makes sense.
  Teuchos::Array< double > work (std::max (1, lwork));
  // Solve the least-squares problem.
  lapack.GELS('N', Vmatrix.numRows(), Vmatrix.numCols(),  RHS.numCols(),
                   Vmatrix.values(),  Vmatrix.numRows(),  RHS.values(),   RHS.numRows(),
                   &work[0], lwork, &info);

  for(int ii=0; ii<basisSize; ii++)  
  {
    samplePCE[ii] = RHS(ii,0); // this line fails if I change the regressionPCE object to type Sacado::PCE::OrthogPoly
  }
}
#endif


#if Xyce_STOKHOS_ENABLE
//-------------------------------------------------------------------------------
// Function      : solveProjectionPCE
// Purpose       : solves for the coefficients of an NISP problem for a single output
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 8/1/2018
//-------------------------------------------------------------------------------
void solveProjectionPCE(
    const Teuchos::RCP<const Stokhos::ProductBasis<int,double> > & basis,
    const Teuchos::RCP<const Stokhos::Quadrature<int,double> > & quad,
    const std::vector<double> & f,
    Sacado::PCE::OrthogPoly<double, Stokhos::StandardStorage<int,double> > & projectionPCE
    )
{
  using Teuchos::Array;
  const int pce_size                        = basis->size();
  const int num_quad_points                 = quad->size();
  const Array<double>& quad_weights         = quad->getQuadWeights();
  const Array< Array<double> >& quad_points = quad->getQuadPoints();
  const Array< Array<double> >& quad_values = quad->getBasisAtQuadPoints();

  TEUCHOS_TEST_FOR_EXCEPTION(num_quad_points != f.size(), std::logic_error,
       "solveProjectionPCE: number of quadrature points"
       << "does not match the output vector size.  f = " << f.size() 
       << ".  num_quad_points = " << num_quad_points << ".");

  // Loop over quadrature points 
  for (int qp_block=0; qp_block<num_quad_points; ++qp_block)
  {
    // Sum results into PCE integral
    for (int pc=0; pc<pce_size; ++pc) 
    {
      const double inv_nrm_sq = 1.0 / basis->norm_squared(pc);

      const double w = quad_weights[qp_block];
      const double psi = quad_values[qp_block][pc];
      projectionPCE.fastAccessCoeff(pc) += inv_nrm_sq * w * f[qp_block] * psi;
    }
  }
}
#endif

//-------------------------------------------------------------------------------
// Function      : applyCovariance
// Purpose       :
//
// This routine applies the covariance matrix and transforms a set of 
// uncorrelated standard normal samples into correlated samples for a given set of
// means and a covariance matrix.
//-------------------------------------------------------------------------------
// X,A,Y are matrices, but are stored in row-major vector form.
// means is a simple vector, one mean per parameter.
//
// A is the covariance matrix
// X, the input, should contain standard normal samples.
// Y, the output will contain correlated samples.
//
// The covariance matrix, A, is a symmetric (m x n) matrix, where m = n = numParams.  
// If it isn't symmetic, then this routine doesn't work.  Being symmetric, 
// it doesn't actually matter if stored as row-major or col-major form.
//
// X and Y are non-symmetric (m x k) matrices, where m = numParams, k = numSamples
// So, each row of X and Y corresponds to a parameter.  Each column corresponds 
// to a different sample.
//
// As X is just uncorrelated random standard normal distributions, the distribution is 
// identical for all parameters, meaning that it doesn't matter if X is row-major or 
// col-major form either.
//
// So, the row-major assumption really only applies to the output Y.
//-------------------------------------------------------------------------------
// Notes on the method:
// This routine ransforms a matrix of standard normals into matrix where each row
// contains multivariate normals with the desired covariance.
//
// Compute L such that dot(transpose(L),L) == covariance.
//
// The matrix products of the rows of x and L will have the desired
// covariance. Sqrt(s)*v where (u,s,v) is the singular value decomposition (SVD) 
// of the covariance matrix is such an L.  So, L = sqrt(S)*V.
// L is the "square root" of the covariance matrix.  Many decompositions give you 
// this, so it doesn't have to be SVD.
//
// The most obvious alternative is to use a Cholesky decomposition.  This might be 
// better b/c Cholesky is often used to determine if a matrix is SPD or not, and 
// as of this writing, this function doesn't check that.  A non-symmetric covariance 
// matrix doesn't make sense.
//
// Another alternative method is to do an Eigenvalue decomposition.  For an 
// SPD matrix, the SVD and Eigenvalue decompositions are the same.  Generally (from 
// what I have read anyway) the SVD decomposition is numerically more stable, 
// but it is slower.  
//
// For now, all my covariance matrices are very small, so
// the choice isn't very important, but that may change later.
//
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 2/14/2018
//-------------------------------------------------------------------------------
void applyCovariance(
    const int numParams,
    const int numSamples,
    const std::vector<double> & X,     // input: uncorrelated samples matrix  (a matrix of standard normals)
    const std::vector<double> & A,     // input: covariance matrix, in row-major storage form
    const std::vector<double> & means, // input: means vector
    std::vector<double> & Y)           // output: correlated samples matrix
{
  Teuchos::LAPACK<int,double> lapack;
  Teuchos::BLAS<int, double> blas;

  const int m = numParams; // using a square matrix so m=n=size.
  const int n = m;
  int lwork=10*m;
  int info=0;
  std::vector<double> U(m*n,0.0), S(m,0.0), VT(m*n,0.0);
  std::vector<double> work(lwork,0.0), rwork(lwork,0.0); 

  std::vector<double> AnonConst = A; // necessary b/c A was passed in as 'const'
  lapack.GESVD('A','A',m,n,&AnonConst[0],m,&S[0],&U[0],m,&VT[0],m,&work[0],lwork,&rwork[0],&info);

  // now that we have the SVD of the covariance matrix, compute L = sqrt(S)*V
  std::vector<double> sqrtS(S.size(),0.0);
  for (int is=0;is<S.size();++is)
  {
    sqrtS[is] = std::sqrt(S[is]);
  }
  std::vector<double> L(A.size(),0.0);

  std::vector<double> sqrtSmatrix(A.size(),0.0);
  for (int row=0;row<m;++row)
  {
    int i = n * row + row; // i = n * row + col, but row=col
    sqrtSmatrix[i] = sqrtS[row]; // put on diagonal
  }

  // compute L.  It can be either U*sqrt(S) or sqrt(S)*V^T
  blas.GEMM( Teuchos::NO_TRANS, Teuchos::NO_TRANS, m, m, m, 1.0,  &U[0], m, &sqrtSmatrix[0], m,  0.0, &L[0],m);

  // BLAS can handle either col-major or row-major storage.  I have been assuming row-major
  // If you use the NO_TRANS option, the referred to object is in col-major format.
  //
  // Since everything is row-major, I have to transpose everything to get right answer.
  // This can be done in more than one way.  The way I've chosen is from this website:
  // https://www.christophlassner.de/using-blas-from-c-with-row-major-data.html
  // If col-major GEMM call is: dgemm("N","N", m, n, k, alpha, A, k, B, n, beta, C, n);
  // then row-major GEMM call is: dgemm("N","N", n, m, k, alpha, B, n, A, k, beta, C, n);
  
  // compute L*X
  blas.GEMM(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 
    numSamples, numParams, numParams, 1.0, &X[0], numSamples, &L[0], numParams, 0.0, &Y[0], numSamples);

  // add in the means.
  for(int col=0;col<numSamples;++col) 
   { for (int row=0;row<numParams;++row) { Y[numSamples * row + col] += means[row]; } }
}


//-------------------------------------------------------------------------------
//
// cephes functions from netlib:  (ndtri, polevl, p1evl)
//
//-------------------------------------------------------------------------------
//							polevl.c
//							p1evl.c
//
//	Evaluate polynomial
//
// int N;
// double x, y, coef[N+1], polevl[];
//
// y = polevl( x, coef, N );
//
//
//
// DESCRIPTION:
//
// Evaluates polynomial of degree N:
//
//                     2          N
// y  =  C  + C x + C x  +...+ C x
//        0    1     2          N
//
// Coefficients are stored in reverse order:
//
// coef[0] = C  , ..., coef[N] = C  .
//            N                   0
//
// The function p1evl() assumes that coef[N] = 1.0 and is
// omitted from the array.  Its calling arguments are
// otherwise the same as polevl().
//
// Cephes Math Library Release 2.1:  December, 1988
// Copyright 1984, 1987, 1988 by Stephen L. Moshier
// Direct inquiries to 30 Frost Street, Cambridge, MA 02140
//
//-------------------------------------------------------------------------------
double polevl( double x, double * coef, int N )
{
  double ans;
  int i;
  double *p;

  p = coef;
  ans = *p++;
  i = N;

  do
    ans = ans * x  +  *p++;
  while( --i );

  return( ans );
}

//-------------------------------------------------------------------------------
//							p1evl()	
//                                          N
// Evaluate polynomial when coefficient of x  is 1.0.
// Otherwise same as polevl.
//-------------------------------------------------------------------------------
double p1evl( double x, double * coef, int N )
{
  double ans;
  double *p;
  int i;

  p = coef;
  ans = x + *p++;
  i = N-1;

  do
    ans = ans * x  + *p++;
  while( --i );

  return( ans );
}

/* sqrt(2pi) */
static double s2pi = 2.50662827463100050242E0;

/* approximation for 0 <= |y - 0.5| <= 3/8 */
static double P0[5] = {
  -5.99633501014107895267E1,
  9.80010754185999661536E1,
  -5.66762857469070293439E1,
  1.39312609387279679503E1,
  -1.23916583867381258016E0,
};

static double Q0[8] = {
  /* 1.00000000000000000000E0, */
  1.95448858338141759834E0,
  4.67627912898881538453E0,
  8.63602421390890590575E1,
  -2.25462687854119370527E2,
  2.00260212380060660359E2,
  -8.20372256168333339912E1,
  1.59056225126211695515E1,
  -1.18331621121330003142E0,
};

/* Approximation for interval z = sqrt(-2 log y ) between 2 and 8
 * i.e., y between exp(-2) = .135 and exp(-32) = 1.27e-14.
 */
static double P1[9] = {
  4.05544892305962419923E0,
  3.15251094599893866154E1,
  5.71628192246421288162E1,
  4.40805073893200834700E1,
  1.46849561928858024014E1,
  2.18663306850790267539E0,
  -1.40256079171354495875E-1,
  -3.50424626827848203418E-2,
  -8.57456785154685413611E-4,
};

static double Q1[8] = {
  /*  1.00000000000000000000E0, */
  1.57799883256466749731E1,
  4.53907635128879210584E1,
  4.13172038254672030440E1,
  1.50425385692907503408E1,
  2.50464946208309415979E0,
  -1.42182922854787788574E-1,
  -3.80806407691578277194E-2,
  -9.33259480895457427372E-4,
};

/* Approximation for interval z = sqrt(-2 log y ) between 8 and 64
 * i.e., y between exp(-32) = 1.27e-14 and exp(-2048) = 3.67e-890.
 */

static double P2[9] = {
  3.23774891776946035970E0,
  6.91522889068984211695E0,
  3.93881025292474443415E0,
  1.33303460815807542389E0,
  2.01485389549179081538E-1,
  1.23716634817820021358E-2,
  3.01581553508235416007E-4,
  2.65806974686737550832E-6,
  6.23974539184983293730E-9,
};

static double Q2[8] = {
  /*  1.00000000000000000000E0, */
  6.02427039364742014255E0,
  3.67983563856160859403E0,
  1.37702099489081330271E0,
  2.16236993594496635890E-1,
  1.34204006088543189037E-2,
  3.28014464682127739104E-4,
  2.89247864745380683936E-6,
  6.79019408009981274425E-9,
};

// replace with something else.... see below
//static double MAXNUM =  1.79769313486231570815E308; 

//-----------------------------------------------------------------------------
// Function      : ndtri
// Purpose       : computes the inverse of a normal distribution
//
// Special Notes : this came from the cephes library, which I got from netlib.
//                 The license is free.
//
// x = ndtri( y );
//
// Returns the argument, x, for which the area under the
// Gaussian probability density function (integrated from
// minus infinity to x) is equal to y.
//
//
// For small arguments 0 < y < exp(-2), the program computes
// z = sqrt( -2.0 * log(y) );  then the approximation is
// x = z - log(z)/z  - (1/z) P(1/z) / Q(1/z).
// There are two rational functions P/Q, one for 0 < y < exp(-32)
// and the other for y up to exp(-2).  For larger arguments,
// w = y - 0.5, and  x/sqrt(2pi) = w + w**3 R(w**2)/S(w**2)).
//
// Cephes Math Library Release 2.8:  June, 2000
// Copyright 1984, 1987, 1989, 2000 by Stephen L. Moshier
//
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 8/15/2018
//-----------------------------------------------------------------------------
double ndtri(double y0)
{
  double x, y, z, y2, x0, x1;
  int code;

  if( y0 <= 0.0 )
  {
    Xyce::Report::DevelFatal0().in("ndtri - ") << "Function out of bounds.";
    return( std::numeric_limits<double>::min() );
  }
  if( y0 >= 1.0 )
  {
    Xyce::Report::DevelFatal0().in("ndtri - ") << "Function out of bounds.";
    return( std::numeric_limits<double>::max() );
  }
  code = 1;
  y = y0;
  if( y > (1.0 - 0.13533528323661269189) ) /* 0.135... = exp(-2) */
  {
    y = 1.0 - y;
    code = 0;
  }

  if( y > 0.13533528323661269189 )
  {
    y = y - 0.5;
    y2 = y * y;
    x = y + y * (y2 * polevl( y2, P0, 4)/p1evl( y2, Q0, 8 ));
    x = x * s2pi; 
    return(x);
  }

  x = std::sqrt( -2.0 * std::log(y) );
  x0 = x - std::log(x)/x;

  z = 1.0/x;
  if( x < 8.0 ) /* y > exp(-32) = 1.2664165549e-14 */
    x1 = z * polevl( z, P1, 8 )/p1evl( z, Q1, 8 );
  else
    x1 = z * polevl( z, P2, 8 )/p1evl( z, Q2, 8 );
  x = x0 - x1;
  if( code != 0 )
    x = -x;
  return( x );
}

//-----------------------------------------------------------------------------
// Function      : inverf
// Purpose       :
// Special Notes : 
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 8/15/2018
//-----------------------------------------------------------------------------
double inverf (double p)
{
  // ndtri will give you the Inverse of Normal distribution function
  // This converts it to the inverse of an error function
  double root_2 = std::sqrt(2);
  double arg    = (p+1)/2.0;
  double retval = ndtri(arg)/root_2;
  return retval;
}

//-----------------------------------------------------------------------------
// Function      : setupNormal
// Purpose       : takes a point from a uniform distribution and converts to a normal distribution
// Special Notes : 
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 9/27/2018
//-----------------------------------------------------------------------------
double setupNormal(const double val, const double mean, const double stddev)
{
  // this part converts the uniform sampling to standard normal distribution via the inverse transform method
  double retVal = ndtri(val);

  // convert from std normal to the given mean and std dev.
  retVal *= stddev;
  retVal += mean;
  return retVal;
}

//-----------------------------------------------------------------------------
// Function      : setupUniform
// Purpose       : converts from uniform distribution from [0,1] to one with a specified min/max
// Special Notes : 
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 9/27/2018
//-----------------------------------------------------------------------------
double setupUniform(const double val, const double min, const double max)
{
  // stretch/compress this uniform distribution to fit the given min and max.
  double dv = fabs(max-min);
  return (dv*val + min);
}

//-----------------------------------------------------------------------------
// Function      : setupMonteCarloSampleValues
// Purpose       :
// Special Notes : 
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 8/15/2018
//-----------------------------------------------------------------------------
void setupMonteCarloSampleValues(
    long theSeed,
    const int numSamples,
    const SweepVector & samplingVector,
    std::vector<double> & X)
{
  int numParams = samplingVector.size();

#if __cplusplus>=201103L
  // allocate the Mersenne Twister algorithm class
  std::mt19937 * mtPtr = new std::mt19937(theSeed);
  std::mt19937 & mt = *mtPtr;

  // sampling related objects:
  std::uniform_real_distribution<double> uniformDistribution(0.0,1.0);
  //std::normal_distribution<double> normalDistribution(1.0,0.0);
#else
  Xyce::Util::RandomNumbers *theRandomNumberGenerator = new Xyce::Util::RandomNumbers(theSeed);
#endif

  X.resize(numSamples*numParams,0.0);

  double val=0.0;

  //bool allDistributionsNormal=true;
  for (int is=0;is<numSamples;++is)
  {
    for (int ip=0;ip<numParams;++ip )
    {
      // calling distribution violates const, but I want most of the sp class to be const in this function
      const SweepParam & sp = samplingVector[ip];

      val=0.0;

#if __cplusplus>=201103L
      if (sp.type == "UNIFORM")
      {
        //allDistributionsNormal=false;
        double prob = uniformDistribution(mt);
        val = setupUniform(prob, sp.startVal, sp.stopVal);
      }
      else if (sp.type == "NORMAL")
      {
        double prob = uniformDistribution(mt);
        val = setupNormal(prob,sp.mean,sp.stdDev);

        while (  (sp.upper_boundGiven && sp.upper_bound < val) ||
                 (sp.lower_boundGiven && sp.lower_bound > val) )
        {
          double prob = uniformDistribution(mt);
          val = setupNormal(prob,sp.mean,sp.stdDev);
        }
      }
      else if (sp.type == "GAMMA")
      {
        //allDistributionsNormal=false;
        std::gamma_distribution<double> gammaDistribution(sp.alpha,sp.beta);
        val = gammaDistribution(mt);

        while (  (sp.upper_boundGiven && sp.upper_bound < val) ||
                 (sp.lower_boundGiven && sp.lower_bound > val) )
        {
          val = gammaDistribution(mt);
        }
      }
#else
      if (sp.type == "UNIFORM")
      {
        //allDistributionsNormal=false;  
        double tmp = theRandomNumberGenerator->uniformRandom();
        double dv = sp.stopVal-sp.startVal;
        val = sp.startVal + dv*tmp;
      }
      else if (sp.type == "NORMAL")
      {
        val = theRandomNumberGenerator->gaussianRandom(sp.mean,sp.stdDev);

        while (  (sp.upper_boundGiven && sp.upper_bound < val) ||
                 (sp.lower_boundGiven && sp.lower_bound > val) )
        {
          val = theRandomNumberGenerator->gaussianRandom(sp.mean,sp.stdDev);
        }
      }
#endif

      X[numSamples * ip + is] = val;
    }
  }

#if __cplusplus>=201103L
  delete mtPtr;
#else
  delete theRandomNumberGenerator;
#endif
}

//-----------------------------------------------------------------------------
// Function      : setupLHSSampleValues
// Purpose       :
// Special Notes : 
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 8/15/2018
//-----------------------------------------------------------------------------
void setupLHSSampleValues(
    long theSeed,
    const int numSamples,
    const SweepVector & samplingVector,
    std::vector<double> & X)
{
#if __cplusplus>=201103L
  // allocate the Mersenne Twister algorithm class
  std::mt19937 * mtPtr = new std::mt19937(theSeed);
  std::mt19937 & mt = *mtPtr;
  std::uniform_int_distribution<> dis(1, numSamples);
  std::uniform_real_distribution<double> distribution(0.0,1.0);
#else
  Xyce::Util::RandomNumbers *theRandomNumberGenerator = new Xyce::Util::RandomNumbers(theSeed);
#endif  
  
  int numParams = samplingVector.size();
  X.resize(numSamples*numParams,0.0);

  for (int ip=0;ip<numParams;++ip )
  {
    const SweepParam & sp = samplingVector[ip];

    for (int is=0;is<numSamples;++is)
    {
      double val=0.0;
      double finalVal=0.0;
      double prob=0.0;
      int bin=0;

#if __cplusplus>=201103L
      val = distribution(mt);
      bin = dis(mt);
#else
      val = theRandomNumberGenerator->uniformRandom(); 
      bin = theRandomNumberGenerator->uniformRandomInt(1, numSamples);
#endif
      prob = (bin-val)/numSamples;

      if (sp.type == "UNIFORM")     
      {
        finalVal = setupUniform(prob, sp.startVal, sp.stopVal);
      }
      else if (sp.type == "NORMAL") 
      {
        finalVal = setupNormal(prob,sp.mean,sp.stdDev);

        while ( (sp.upper_boundGiven && sp.upper_bound < finalVal) ||
                (sp.lower_boundGiven && sp.lower_bound > finalVal) )
        {
#if __cplusplus>=201103L
          val = distribution(mt);
          bin = dis(mt);
#else
          val = theRandomNumberGenerator->uniformRandom(); 
          bin = theRandomNumberGenerator->uniformRandomInt(1, numSamples);
#endif
          prob = (bin-val)/numSamples;
          finalVal = setupNormal(prob,sp.mean,sp.stdDev);
        }
      }
      else
      {
        Xyce::Report::DevelFatal0().in(" setupLHSSampleValues - ") << 
         sp.type << "is an unsupported distribution for LHS.";
      }

      X[numSamples * ip + is] = finalVal;
    }
  }

#if __cplusplus>=201103L
  delete mtPtr;
#else
  delete theRandomNumberGenerator;
#endif
}

//-----------------------------------------------------------------------------
// Function      : setupMonteCarloStdNormals
// Purpose       :
// Special Notes : 
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 9/27/2018
//-----------------------------------------------------------------------------
void setupMonteCarloStdNormals(long theSeed, std::vector<double> & X)
{
#if __cplusplus>=201103L
  // allocate the Mersenne Twister algorithm class
  std::mt19937 * mtPtr = new std::mt19937(theSeed);
  std::mt19937 & mt = *mtPtr;

  // setup X with std-normals: 
  std::normal_distribution<double> stdNormalDistribution(0.0,1.0);
  for (int i=0;i<X.size();++i)
  {
    X[i] = stdNormalDistribution(mt);
  }
  delete mtPtr;
#else
  Xyce::Util::RandomNumbers *theRandomNumberGenerator = new Xyce::Util::RandomNumbers(theSeed);
  for (int i=0;i<X.size();++i)
  {
    X[i] = theRandomNumberGenerator->gaussianRandom(0.0,1.0);
  }
  delete theRandomNumberGenerator;
#endif
}

//-----------------------------------------------------------------------------
// Function      : setupLHSStdNormals
// Purpose       :
// Special Notes : 
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 9/27/2018
//-----------------------------------------------------------------------------
void setupLHSStdNormals(long theSeed, const int numSamples, std::vector<double> & X)
{
  double val=0.0;
  double prob=0.0;
  int bin=0;

#if __cplusplus>=201103L
  // allocate the Mersenne Twister algorithm class
  std::mt19937 * mtPtr = new std::mt19937(theSeed);
  std::mt19937 & mt = *mtPtr;
  std::uniform_int_distribution<> dis(1, numSamples);
  std::uniform_real_distribution<double> distribution(0.0,1.0);

  // setup X with std-normals: 
  for (int i=0;i<X.size();++i)
  {
    val = distribution(mt);
    bin = dis(mt);
    prob = (bin-val)/numSamples;
    X[i] = setupNormal(prob,0.0,1.0);
  }
  delete mtPtr;
#else
  Xyce::Util::RandomNumbers *theRandomNumberGenerator = new Xyce::Util::RandomNumbers(theSeed);
  for (int i=0;i<X.size();++i)
  {
    val = theRandomNumberGenerator->uniformRandom(); 
    bin = theRandomNumberGenerator->uniformRandomInt(1, numSamples);
    prob = (bin-val)/numSamples;
    X[i] = setupNormal(prob,0.0,1.0);
  }
  delete theRandomNumberGenerator;
#endif
}

//-----------------------------------------------------------------------------
// Function      : setupSampleValues
// Purpose       :
// Special Notes : 
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 2/19/2018
//-----------------------------------------------------------------------------
void setupSampleValues(
    long theSeed,
    const UQ::SampleType sampleType,
    const int numSamples,
    const SweepVector & samplingVector,  // STL vector of "sweep param" class objects, one for each parameter
    const std::vector<double> & covMatrix,// input: covariance matrix, in row-major storage form
    const std::vector<double> & meanVec,  // input: means vector
    std::vector<double> & X,              // output: uncorrelated samples matrix 
    std::vector<double> & Y)              // output: correlated samples matrix
{
  int numParams = samplingVector.size();
  X.resize(numSamples*numParams,0.0);
  Y.resize(numSamples*numParams,0.0);

  // if applying covariance, then all the distributions need to be normal
  // distributions, and I don't have a way to apply upper and lower bounds
  // to the distribution.
  bool covMatrixGiven = !(covMatrix.empty());
  if ( covMatrixGiven )
  {
    const int m = numParams;
    const int n = m;

    if (sampleType == UQ::MC)
    {
      setupMonteCarloStdNormals(theSeed,X);
    }
    else if (sampleType == UQ::LHS)
    {
      setupLHSStdNormals(theSeed,numSamples,X);
    }

    applyCovariance(numParams, numSamples, X, covMatrix, meanVec, Y);
  }
  else
  {
    if (sampleType == UQ::MC)
    {
      setupMonteCarloSampleValues(theSeed, numSamples, samplingVector, X);
    }
    else if (sampleType == UQ::LHS)
    {
      setupLHSSampleValues(theSeed, numSamples, samplingVector, X);
    }

    // no covariance to apply, so Y=X
    Y = X;
  }
}


//-----------------------------------------------------------------------------
// Function      : unScaleSampleValues
// Purpose       :
//
// set up the lower-case "x" vector.  Capital "X" (or "Y") contains the 
// actual parameter values used by the device models.
//
// Lower case "x" contains the unscaled values used by PCE.
// So, for Hermite polynomials (normal dist), they must be converted to standard normals.
// and for Legendre polynomials (uniform dist), they must be converted to domain and range over (-1,+1)
//
// This x-vector is needed by the regression PCE algorithm, to set up the Vandermode matrix.
// It is also needed by PCE algorithms in order to "resample" a solved PCE expansion.   
// The random samples that are inputs to the PCE expansion need to be in this unscaled form.
//
// The input, Y_, is a compressed row storage 1D vector, while the output, x, is a 2D array.
// Y_ is generated by the function setupSampleValues.  That function requires compressed row storage, 
// but only because the function applyCovariance requires it (and that function isn't always used).
//
// The 2D storage form of x is more convenient, so I generally prefer to use it unless I have no choice.
//
// So, for that reason, the functions in this file a inconsistent in this way.
//
// Also this function will not work correctly if a covariance matrix was applied to create 
// the original scaled values. FIX THIS!
//
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 2/8/2019
//-----------------------------------------------------------------------------
void unScaleSampleValues(
    const int numSamples,
    const SweepVector & samplingVector,
    const std::vector<double> & covMatrix,// input: covariance matrix, in row-major storage form
    const std::vector<double> & meanVec,  // input: means vector
    const std::vector<double> & Y, 
    std::vector< std::vector<double> > & x
    )
{
  const int numParams = samplingVector.size();
  x.resize(numParams);
  for (int i=0;i<numParams;i++) { x[i].resize(numSamples ,0.0);}

  bool covMatrixGiven = !(covMatrix.empty());

  for(int col=0;col<numSamples;++col)
  {
    for (int row=0;row<numParams;++row)
    { 
      const SweepParam & sp = samplingVector[row];

      if ( covMatrixGiven )  
      {
        // only works if everything is normal, but this should have been tested earlier
        double sigma = std::sqrt(fabs(covMatrix[row*numParams+row]));
        double mu = meanVec[row]; 
        x[row][col] = (1.0/sigma) * (Y[numSamples * row + col] - mu);
      }
      else
      {
        if (sp.type == "NORMAL")
        {
          double mu=sp.mean;
          double sigma=sp.stdDev;
          x[row][col] = (1.0/sigma) * (Y[numSamples * row + col] - mu);
        }
        else if (sp.type == "UNIFORM")
        {
          double max = sp.stopVal;
          double min = sp.startVal;
          double delta = fabs(max-min);
          double xSlope = (Y[numSamples * row + col] - min)/delta;
          x[row][col]  = 2.0*xSlope - 1.0;
        }
      }
    }
  }
}

#if Xyce_STOKHOS_ENABLE
//-----------------------------------------------------------------------------
// Function      : setupPCEQuadPoints
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 9/4/2018
//-----------------------------------------------------------------------------
void setupPCEQuadPoints (
    const Teuchos::RCP<const Stokhos::ProductBasis<int,double> > & basis,
    const Teuchos::RCP<const Stokhos::Quadrature<int,double> > & quad,
    const Teuchos::RCP<Stokhos::QuadOrthogPolyExpansion<int,double> > & expn,
    const SweepVector & samplingVector,  // STL vector of "sweep param" class objects, one for each parameter
    const std::vector<double> & covMatrix,// input: covariance matrix, in row-major storage form
    const std::vector<double> & meanVec,  // input: means vector
    std::vector<double> & X,              // output: uncorrelated points matrix
    std::vector<double> & Y)              // output: correlated points matrix
{
  using Teuchos::Array;
  using Teuchos::RCP;
  using Teuchos::rcp;

  int numParams = samplingVector.size();

  bool covMatrixGiven = !(covMatrix.empty());
  if ( covMatrixGiven )
  //if ( covMatrixGiven && allDistributionsNormal )
  {
    // not supported here...
    Xyce::dout() << "Error: cov matrix can't work yet with quadrature PCE" << std::endl;
    exit(0);
  }
  else
  {
    ////////////////////////////////////////////////////////////////////////////////
    // now compute it projection PCE (NISP)
    // Typename of Polynomial Chaos scalar type
    typedef Stokhos::StandardStorage<int,double> pce_storage_type;
    typedef Sacado::PCE::OrthogPoly<double, pce_storage_type> pce_type;

    //
    // Compute PCE expansion of function using NISP 
    // Extract quadrature data
    const int pce_size                        = basis->size();
    const int num_quad_points                 = quad->size();
    const Array<double>& quad_weights         = quad->getQuadWeights();
    const Array< Array<double> >& quad_points = quad->getQuadPoints();
    const Array< Array<double> >& quad_values = quad->getBasisAtQuadPoints();
    const int numSamples = num_quad_points;
    X.resize(numSamples*numParams,0.0);
    Y.resize(numSamples*numParams,0.0);

    // Loop over quadrature points 
    for (int qp_block=0; qp_block<num_quad_points; ++qp_block)
    {
      for (int ip=0;ip<numParams;++ip )
      {
        const SweepParam & sp = samplingVector[ip];
        double val = 0.0;
        pce_type param_pce(expn);

        if (sp.type == "UNIFORM")
        {
          param_pce.term(ip,0) = (sp.stopVal + sp.startVal)/2.0; // mean
          param_pce.term(ip,1) = (sp.stopVal - sp.startVal)/2.0; // slope
          val = param_pce.evaluate(quad_points[qp_block], quad_values[qp_block]);
        }
        else if (sp.type == "NORMAL") 
        {
          param_pce.term(0,0) = sp.mean;
          param_pce.term(ip,1) = sp.stdDev;
          val = param_pce.evaluate(quad_points[qp_block], quad_values[qp_block]);
        }

        int is=qp_block;
        X[numSamples * ip + is] = val;
      }
    }

    Y = X;
  }

}
#endif


//-----------------------------------------------------------------------------
// Function      : checkParameterList
// Purpose       : checks the list of params to ensure they all exist
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 2/19/2018
//-----------------------------------------------------------------------------
void checkParameterList(
    Parallel::Machine comm, 
    Loader::Loader &loader, 
    std::vector<SweepParam>::iterator begin, 
    std::vector<SweepParam>::iterator end)
{
  // loop over the param containers, and check that all the params exist.
  // The device package will complain and exit if it can't find the param.
  for (std::vector<SweepParam>::iterator it = begin; it != end; ++it)
  {
    SweepParam &sweep_param = (*it);
    loader.getParamAndReduce(comm, sweep_param.name);
  }
}

//-----------------------------------------------------------------------------
// Function      : updateSamplingParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, 
// Creation Date : 8/14/2017
//-----------------------------------------------------------------------------
bool updateSamplingParams(
    Loader::Loader &loader, 
    int sample, 
    std::vector<SweepParam>::iterator begin, 
    std::vector<SweepParam>::iterator end, 
    const std::vector<double> & Y,
    int numSamples,
    bool overrideOriginal)
{
  bool reset = false;

  // set parameter(s)
  int ip=0;
  for (std::vector<SweepParam>::iterator it = begin; it != end; ++it,++ip)
  {
    (*it).currentVal = Y[numSamples * ip + sample];
    std::string setParamName;
    Xyce::Nonlinear::getSetParamName( (*it).name, setParamName);
    loader.setParam(setParamName, (*it).currentVal, overrideOriginal);
  }

  return reset;
}

//-----------------------------------------------------------------------------
// Function      : getTheSeed
// Purpose       : Deal with the random number seed.  
// Special Notes :
//
// Usually the user will not set this.  However, for reproducible testing, it can be
// necessary for the user to set the seed.  The user can set the seed either from
// the command line or from the netlist.  Command line overrides the netlist.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, 
// Creation Date : 8/31/2018
//-----------------------------------------------------------------------------
long getTheSeed(
    Parallel::Machine comm, 
    const Xyce::IO::CmdParse & commandLine, int userSeed, bool userSeedGiven)
{
#if __cplusplus>=201103L
  std::random_device rd;
#endif

  long theSeed;
  if (commandLine.argExists("-randseed"))
  {
    long cmdSeed;
    std::stringstream iss(commandLine.getArgumentValue("-randseed"));
    iss >> cmdSeed;
    theSeed = cmdSeed;
  }
  else
  {
    if(userSeedGiven)
    {
      theSeed = userSeed;
    }
    else
    {
      theSeed = 0;
      if (Parallel::is_parallel_run(comm)) 
      {
        N_ERH_ErrorMgr::safeBarrier(comm);

        if (Parallel::rank(comm) == 0) 
        {
#if __cplusplus>=201103L
          theSeed = rd();
#else
          theSeed=time(NULL);
#endif
          Util::Marshal mout;
          mout << theSeed;
          Parallel::Broadcast(comm, mout, 0);
        }
        else 
        {
          Util::Marshal min("");
          Parallel::Broadcast(comm, min, 0);
          min >> theSeed;
        }
      }
      else
      {
#if __cplusplus>=201103L
        theSeed = rd();
#else
        theSeed=time(NULL);
#endif
      }
    }
  }

#if __cplusplus>=201103L
  Xyce::lout() << "Seeding random number generator with " << theSeed << std::endl;
#endif

  return theSeed;
}

#if 0
//-----------------------------------------------------------------------------
// print a histrogram.  doesn't work yet 
//-----------------------------------------------------------------------------
void histrogram( 
    //std::ofstream & fout, 
    std::ostream & fout, 
    const std::string & name, 
    const std::vector<double> function )
{
  int numSamples = function.size();

  std::map<int, int> hist{};
  for(int is=0; is<numSamples; ++is)
  {
    ++hist[std::round(function[is])];
  }
  //int dS = numSamples/50;
  int dS = 10;

  if (dS<=0) {dS = 1;}

  for(auto p : hist)
  {
    fout << std::setw(2)
         << p.first << ' ' << std::string(p.second/dS, '*') << '\n';
  }
  fout << std::endl;
}
#endif

} // namespace UQ
} // namespace Analysis
} // namespace Xyce

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


#ifndef Xyce_N_ANP_UQSupport_h
#define Xyce_N_ANP_UQSupport_h

#include <N_PDS_fwd.h>
#include <N_IO_fwd.h>
#include <N_ANP_UQ_fwd.h>
#include <N_ANP_SweepParam.h>
#include <N_ANP_SweepParamFreeFunctions.h>

#include <N_UTL_ExpressionData.h>
#include <N_UTL_MachDepParams.h>

#ifdef Xyce_STOKHOS_ENABLE
#include <Stokhos_Sacado.hpp>
#include <Sacado_No_Kokkos.hpp>
#endif

#include <vector>

namespace Xyce {
namespace Analysis {
namespace UQ {

//-----------------------------------------------------------------------------
struct statisticalMoments
{
  statisticalMoments() :
    mean(0.0), stddev(0.0), variance(0.0), skew(0.0), kurtosis(0.0),
    max(-Xyce::Util::MachineDependentParams::DoubleMax()),
    min(+Xyce::Util::MachineDependentParams::DoubleMax())
  {};

  double mean;
  double stddev;
  double variance;
  double skew;
  double kurtosis;
  double max;
  double min;
};

//-----------------------------------------------------------------------------
void computeStats (const std::vector<double> & values, statisticalMoments & sm);

//-----------------------------------------------------------------------------
inline std::ostream & operator<<(std::ostream & os, const UQ::statisticalMoments & sm) 
{
  os << "mean     = " << sm.mean <<std::endl;
  os << "stddev   = " << sm.stddev <<std::endl;
  os << "variance = " << sm.variance <<std::endl;
  os << "skew     = " << sm.skew <<std::endl;
  os << "kurtosis = " << sm.kurtosis <<std::endl;
  os << "max      = " << sm.max <<std::endl;
  os << "min      = " << sm.min <<std::endl;

  return os;
}

//-----------------------------------------------------------------------------
// Class         : output function data
//
// Purpose       : This class contains information associated with 
//                 each requested UQ output.
//
// Creator       : Eric Keiter, SNL
// Creation Date : 8/17/2015
//-----------------------------------------------------------------------------
class outputFunctionData
{
public:  
  outputFunctionData() : expDataPtr(0), 
  measureResponseFound(false) {};
  ~outputFunctionData() { if (expDataPtr!=0) { delete expDataPtr; } }

public:
  std::vector<int> globalParamVariableStencil;
  std::string outFuncString;
  Util::ExpressionData * expDataPtr;
  bool measureResponseFound;
  statisticalMoments sm;
  std::vector<double> sampleOutputs;
#ifdef Xyce_STOKHOS_ENABLE
  Stokhos::OrthogPolyApprox<int,double> regressionPCE;
  Sacado::PCE::OrthogPoly<double, Stokhos::StandardStorage<int,double> > projectionPCE;
#endif

  void completeStatistics ()
  {
    computeStats (sampleOutputs, sm);
  };

  void output(std::ostream & os, std::string type="") const
  {
    os << std::endl;
    os << type << " sampling mean of " << outFuncString << " = " << sm.mean <<std::endl;
    os << type << " sampling stddev of " << outFuncString << " = " << sm.stddev <<std::endl;
    os << type << " sampling variance of " << outFuncString << " = " << sm.variance <<std::endl;
    os << type << " sampling skew of " << outFuncString << " = " << sm.skew <<std::endl;
    os << type << " sampling kurtosis of "<< outFuncString << " = " << sm.kurtosis <<std::endl;
    os << type << " sampling max of "<< outFuncString << " = " << sm.max <<std::endl;
    os << type << " sampling min of "<< outFuncString << " = " << sm.min <<std::endl;
  }
};

//-----------------------------------------------------------------------------
// Function      : output function data operator<<
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date :
//-----------------------------------------------------------------------------
inline std::ostream & operator<<(std::ostream & os, const outputFunctionData & ofd)
{
  ofd.output(os);
  return os;
}

#ifdef Xyce_STOKHOS_ENABLE
void solveRegressionPCE(
    const int numParams, // dimension (number of params)
    const std::vector< std::vector<double> > & x,
    const std::vector<double> & f,
    Stokhos::OrthogPolyApprox<int,double> & samplePCE
    );
#endif

#ifdef Xyce_STOKHOS_ENABLE
void solveProjectionPCE(
    const Teuchos::RCP<const Stokhos::ProductBasis<int,double> > & basis,
    const Teuchos::RCP<const Stokhos::Quadrature<int,double> > & quad,
    const std::vector<double> & f,
    Sacado::PCE::OrthogPoly<double, Stokhos::StandardStorage<int,double> > & projectionPCE
    );
#endif

void setupMonteCarloStdNormals(long theSeed, std::vector<double> & X);
void setupLHSStdNormals(long theSeed, const int numSamples, std::vector<double> & X);

void applyCovariance(
    const int numParams,
    const int numSamples,
    const std::vector<double> & X,     // input: uncorrelated samples matrix 
    const std::vector<double> & A,     // input: covariance matrix, in row-major storage form
    const std::vector<double> & means, // input: means vector
    std::vector<double> & Y);          // output: correlated samples matrix

double polevl(double x, double * coef, int N);
double p1evl(double x, double * coef, int N);
double ndtri(double y0);

double setupNormal(const double val, const double mean, const double stddev);
double setupUniform(const double val, const double min, const double max);

void setupMonteCarloSampleValues(
    long theSeed,
    const int numSamples,
    const SweepVector & samplingVector,
    std::vector<double> & X);

void setupLHSSampleValues(
    long theSeed,
    const int numSamples,
    const SweepVector & samplingVector,
    std::vector<double> & X);

void setupSampleValues(
    long theSeed,
    const UQ::SampleType sampleType,
    const int numSamples,
    const SweepVector & samplingVector,
    const std::vector<double> & covMatrix,// input: covariance matrix, in row-major storage form
    const std::vector<double> & meanVec,  // input: means vector
    std::vector<double> & X,              // output: uncorrelated samples matrix 
    std::vector<double> & Y);              // output: correlated samples matrix

void unScaleSampleValues(
    const int numSamples,
    const SweepVector & samplingVector,
    const std::vector<double> & covMatrix,
    const std::vector<double> & meanVec,
    const std::vector<double> & Y, 
    std::vector< std::vector<double> > & x
    );

#ifdef Xyce_STOKHOS_ENABLE
void setupPCEQuadPoints (
    const Teuchos::RCP<const Stokhos::ProductBasis<int,double> > & basis,
    const Teuchos::RCP<const Stokhos::Quadrature<int,double> > & quad,
    const Teuchos::RCP<Stokhos::QuadOrthogPolyExpansion<int,double> > & expn,
    const SweepVector & samplingVector,
    const std::vector<double> & covMatrix,// input: covariance matrix, in row-major storage form
    const std::vector<double> & meanVec,  // input: means vector
    std::vector<double> & X,              // output: uncorrelated points matrix 
    std::vector<double> & Y);              // output: correlated points matrix
#endif

void checkParameterList(
    Parallel::Machine comm, 
    Loader::Loader &loader, 
    std::vector<SweepParam>::iterator begin, 
    std::vector<SweepParam>::iterator end);

bool updateSamplingParams(
    Loader::Loader &loader, 
    int sample, 
    std::vector<SweepParam>::iterator begin, 
    std::vector<SweepParam>::iterator end, 
    const std::vector<double> & Y,
    int numSamples,
    bool overrideOriginal);

bool updateExpressionSamplingTerms2(
    Loader::Loader &loader, 
    int sample, 
    std::vector<SweepParam> & samplingVec,
    const std::vector<double> & Y,
    int numSamples,
    bool overrideOriginal);

long getTheSeed(
    Parallel::Machine comm, 
    const Xyce::IO::CmdParse & commandLine, 
    int userSeed, 
    bool userSeedGiven,
    bool output=true);

#ifdef Xyce_STOKHOS_ENABLE
//-----------------------------------------------------------------------------
// Function      : sampleApproximationPCE
//
// Purpose       : Do sampling for a list of uncertain inputs and outputs using 
//                 a pre-computed PCE approximation.  Once this is done, 
//                 compute sampling-based statistics such as mean, variance, 
//                 kurtosis, etc.
//
//                 The number of sample points is generally going to be a lot 
//                 larger than the number of evaluation points used to solve for 
//                 the PCE coefficients.   So, a new set of sample points for the
//                 list of uncertain inputs must be generated first, then 
//                 unscaled, and then used as input to the PCE approximation.
// 
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 2/10/2019
//-------------------------------------------------------------------------------
template <typename ScalarT>
void sampleApproximationPCE(
    long theSeed,
    const UQ::SampleType sampleType,
    const SweepVector & samplingVector,
    const std::vector<double> & covMatrix,
    const std::vector<double> & meanVec,
    const int numSamples,
    const int numParams,
    const std::vector< ScalarT > & samplePCEvec,
    std::vector <std::vector<double> > & fvec,
    std::vector <UQ::statisticalMoments> & statVec)
{
  // 0. check size of samplePCEvec and fvec (should both = number of outputs)
  // 1. allocate sampling objects, X_, Y_, and x.
  // 2. setup sample values in X_ and Y_ using the setupSampleValues function
  // 3. unScale it to create x.
  // 4. loop over x, and create points, which are fed into the samplePCE object.
  // 5. place the results of those evaluations over x into the f object.
  // 6. compute statistics

  // 0.
  if (samplePCEvec.size() == 0) return;

  TEUCHOS_TEST_FOR_EXCEPTION(samplePCEvec.size() != fvec.size(), std::logic_error,
       "sampleApproximationPCE: size of f vector "
       << "does not match the PCE vector size.  f.size = " << fvec.size() 
       << ".  PCE.size = " << samplePCEvec.size() << ".");

  // 1.
  int numOutputs=fvec.size();
  for (int ii=0;ii<numOutputs;++ii) { fvec[ii].resize(numSamples,0.0); }

  std::vector<double> X(numSamples*numParams,0.0);
  std::vector<double> Y(numSamples*numParams,0.0);
  std::vector< std::vector<double> > x(numParams);
  for (int i=0;i<numParams;i++) { x[i].resize(numSamples ,0.0);}

  // 2.
  setupSampleValues(theSeed, sampleType, numSamples, samplingVector, covMatrix, meanVec, X, Y);

  // 3.
  unScaleSampleValues(numSamples, samplingVector, covMatrix, meanVec, Y, x);

  // 4, 5.
  for(int iout=0;iout<numOutputs;++iout)
  {
    std::vector<double> & f1d = fvec[iout];
    statisticalMoments & sm = statVec[iout];

    for(int col=0;col<numSamples;++col)
    { 
      Teuchos::Array<double> point( numParams );
      for (int row=0;row<numParams;++row)
      { 
        point[row] = x[row][col];
      }
   
      f1d[col] = samplePCEvec[iout].evaluate(point);
    } // resample loop

    computeStats (f1d, sm);

  }//outputs loop
}

//-----------------------------------------------------------------------------
// Function      : evaluateApproximationPCE
//
// Purpose       : This is a similar function to sampleApproximationPCE, except 
//                 that instead of evaluating the PCE approximation using 
//                 randomly sampled inputs, it instead evaluates it with a 
//                 specific set of parameter values.  
//
//                 It should basically do the same steps 3,4,5 in  
//                 function sampleApproximationPCE.
//
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 8/26/2019
//-------------------------------------------------------------------------------
template <typename ScalarT>
void evaluateApproximationPCE(
    const SweepVector & samplingVector,
    const std::vector<double> & Y,
    const int numSamples,
    const std::vector< ScalarT > & samplePCEvec,
    std::vector <std::vector<double> > & fvec)
{
  // 3. unScale it to create x.
  // 4. loop over x, and create points, which are fed into the samplePCE object.
  // 5. place the results of those evaluations over x into the f object.

  if (samplePCEvec.size() == 0) return;

  TEUCHOS_TEST_FOR_EXCEPTION(samplePCEvec.size() != fvec.size(), std::logic_error,
       "sampleApproximationPCE: size of f vector "
       << "does not match the PCE vector size.  f.size = " << fvec.size() 
       << ".  PCE.size = " << samplePCEvec.size() << ".");

  int numOutputs=fvec.size();
  for (int ii=0;ii<numOutputs;++ii) { fvec[ii].resize(numSamples,0.0); }

  int numParams = samplingVector.size();
  std::vector< std::vector<double> > x(numParams);
  for (int i=0;i<numParams;i++) { x[i].resize(numSamples ,0.0);}

  const std::vector<double> covMatrix; // empty.  This just here as a placeholder
  const std::vector<double> meanVec; // empty.  This just here as a placeholder
  unScaleSampleValues(numSamples, samplingVector, covMatrix, meanVec, Y, x);

  for(int iout=0;iout<numOutputs;++iout)
  {
    std::vector<double> & f1d = fvec[iout];

    for(int col=0;col<numSamples;++col)
    { 
      Teuchos::Array<double> point( numParams );
      for (int row=0;row<numParams;++row)
      { 
        point[row] = x[row][col];
      }
   
      f1d[col] = samplePCEvec[iout].evaluate(point);
    } // resample loop
  }//outputs loop
}
#endif

} // namespace UQ
} // namespace Analysis
} // namespace Xyce
#endif // Xyce_N_ANP_UQSupport_h

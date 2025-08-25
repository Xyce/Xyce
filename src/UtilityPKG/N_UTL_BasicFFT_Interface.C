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
// Purpose        : This file contains specializations for the BasicFFT interface
//                  for various vector types.
//
// Special Notes  : 
//
// Creator        : Richard Schiek
//
// Creation Date  : 1/21/25
//
//
//
//
//-------------------------------------------------------------------------
// ---------- Standard Includes ----------

#include <Xyce_config.h>

#include <N_UTL_BasicFFT_Interface.hpp>

// ----------   Other Includes   ----------

#include <iostream>
#include <vector>
#include <N_UTL_Math.h>
#include <complex>

// ---------- Structure definitions ----------

// forward declaration of internal routine.
void calcCTForwardFFTInt( std::vector<double> & inVec, unsigned int numPoints);
void calcCTInverseFFTInt( std::vector<double> & inVec, unsigned int numPoints);

void calcCTForwardFFTComplex( std::vector<std::complex<double> > & inVec);
  
template<>
void N_UTL_BasicFFT_Interface<std::vector<double> >::calculateDFT()
{
  if(nextLargestPowerOf2_ == signalLength_)
  {
    /*
    // uses the complex data type implementation.
    std::vector<std::complex<double> > tempVec(signalLength_);
    for( auto j=0; j<signalLength_; j++)
    {
      tempVec[j] = std::complex<double>((*dftInData_)[j], 0.0 );
    }
    calcCTForwardFFTComplex( tempVec );
    for(auto j=0; j<((*dftOutData_).size()/2); j++)
    {
      (*dftOutData_)[2*j]=tempVec[j].real();
      (*dftOutData_)[2*j+1]=tempVec[j].imag();
    }
    */
    
    std::vector<double> tempVec(2*(*dftInData_).size());
    for( auto j=0; j<(*dftInData_).size(); j++)
    {
      tempVec[2*j]=(*dftInData_)[j];
      tempVec[2*j+1]=0.0;
    }
    calcCTForwardFFTInt( tempVec, tempVec.size() );
    for(auto j=0; j<(*dftOutData_).size(); j++)
    {
      (*dftOutData_)[j]=tempVec[j];
    }
  }
  else
  {
    for( auto k=0; k<(1+signalLength_/2); k++)
    {
      (*dftOutData_)[2*k] = 0.0;
      (*dftOutData_)[2*k+1] = 0.0;
      for( auto j=0; j<signalLength_; j++)
      {
        (*dftOutData_)[2*k] += std::cos(-2*M_PI*j*k/signalLength_)*(*dftInData_)[j];
        (*dftOutData_)[2*k+1] += std::sin(-2*M_PI*j*k/signalLength_)*(*dftInData_)[j];
      }
    }
  }
}

// Calculate IFT with the vectors that have been registered.
template<>
void N_UTL_BasicFFT_Interface<std::vector<double> >::calculateIFT()
{ 
  if(nextLargestPowerOf2_ == signalLength_)
  {
    /*
    // Slower, but simple implementation  
    for( auto k=0; k<signalLength_; k++)
    {
      (*iftOutData_)[k] = 0.0;
      for( auto j=0; j<((signalLength_/2)); j++)
      {
        (*iftOutData_)[k] += (std::cos(2*M_PI*j*k/signalLength_))*(*iftInData_)[2*j] - (std::sin(2*M_PI*j*k/signalLength_))*(*iftInData_)[2*j+1];
      }
      
      for( auto j=(signalLength_/2), jc=(signalLength_/2); j>0; j--, jc++)
      {
        (*iftOutData_)[k] += (std::cos(2*M_PI*jc*k/signalLength_))*(*iftInData_)[2*j] + (std::sin(2*M_PI*jc*k/signalLength_))*(*iftInData_)[2*j+1];
      }
      (*iftOutData_)[k] = (*iftOutData_)[k] / signalLength_;
    }
    */
    
    
    // copy input data to working vector with complex conjugate in upper half of working vector
    std::vector<double> tempVec(2*(*iftInData_).size()-4, 0.0);
    for( auto j=0; j<(((*iftInData_).size()/2)); j++)
    {
      tempVec[2*j]=(*iftInData_)[2*j];
      tempVec[2*j+1]=(*iftInData_)[2*j+1];
      if( j < (((*iftInData_).size()/2)-2))
      {
        tempVec[2*j+(*iftInData_).size()]=(*iftInData_)[(*iftInData_).size()-2*j-4];
        tempVec[2*j+(*iftInData_).size()+1]=-(*iftInData_)[(*iftInData_).size()-2*j-3];
      }
    }
    calcCTInverseFFTInt( tempVec, tempVec.size() );
    for(auto j=0; j<(*iftOutData_).size(); j++)
    {
      (*iftOutData_)[j]=tempVec[2*j]/(*iftOutData_).size();
    }
  }
  else
  {
    for( auto k=0; k<signalLength_; k++)
    {
      (*iftOutData_)[k] = 0.0;
      for( auto j=0; j<((signalLength_/2)+1); j++)
      {
        (*iftOutData_)[k] += (std::cos(2*M_PI*j*k/signalLength_))*(*iftInData_)[2*j] - (std::sin(2*M_PI*j*k/signalLength_))*(*iftInData_)[2*j+1];
      }
      
      for( auto j=(signalLength_/2), jc=(signalLength_/2)+1; j>0; j--, jc++)
      {
        (*iftOutData_)[k] += (std::cos(2*M_PI*jc*k/signalLength_))*(*iftInData_)[2*j] + (std::sin(2*M_PI*jc*k/signalLength_))*(*iftInData_)[2*j+1];
      }
      (*iftOutData_)[k] = (*iftOutData_)[k] / signalLength_;
    }
  }
}

//
// a Cooley-Tukey forward FFT algorithm of O(n log(n)) number of operations for a
// power of 2 number of sample points.  Also assumes that inVec is real and outVec is
// complex with staggared real, imaginary components.


void calcCTForwardFFTInt( std::vector<double> & inVec, unsigned int numPoints)
{
  auto numVals = inVec.size();  
  if( numVals <= 2)
  {
    // trivial case just return 
  }
  else
  {
    // split the input vector 
    std::vector<double> evenVec(numVals/2);
    std::vector<double> oddVec(numVals/2);
    for( auto j=0; j< (numVals/4); j++)
    {
      evenVec[2*j] = inVec[4*j];
      evenVec[2*j+1] = inVec[4*j+1];
      oddVec[2*j] = inVec[4*j+2];
      oddVec[2*j+1] = inVec[4*j+3];
    }
    calcCTForwardFFTInt( evenVec, (numVals/2));
    calcCTForwardFFTInt( oddVec, (numVals/2));
    // now combine even and odd elements.
    for( auto k = 0; k<(numVals/4); k++)
    {
      auto val1real = evenVec[2*k];
      auto val1imag = evenVec[2*k+1];
      auto val2real = std::cos(-4*M_PI*k/numVals)*oddVec[2*k] - std::sin(-4*M_PI*k/numVals)*oddVec[2*k+1];
      auto val2imag = std::sin(-4*M_PI*k/numVals)*oddVec[2*k] + std::cos(-4*M_PI*k/numVals)*oddVec[2*k+1];
      inVec[2*k] = val1real + val2real;
      inVec[2*k+1] = val1imag + val2imag;
      inVec[2*k+(numVals/2)] = val1real - val2real;
      inVec[2*k+(numVals/2)+1]=  val1imag - val2imag;
    }
  }
}


void calcCTInverseFFTInt( std::vector<double> & inVec, unsigned int numPoints)
{
  auto numVals = inVec.size();
  if( numVals <= 2)
  {
    // trivial case of size 1
  }
  else
  {
    // split the input vector 
    std::vector<double> evenVec(numVals/2);
    std::vector<double> oddVec(numVals/2);
    for( auto j=0; j< (numVals/4); j++)
    {
      evenVec[2*j] = inVec[4*j];
      evenVec[2*j+1] = inVec[4*j+1];
      oddVec[2*j] = inVec[4*j+2];
      oddVec[2*j+1] = inVec[4*j+3];
    }
    calcCTInverseFFTInt( evenVec, (numVals/2));
    calcCTInverseFFTInt( oddVec, (numVals/2));
    // now combine even and odd elements.
    for( auto k = 0; k<(numVals/4); k++)
    {
      auto val1real = evenVec[2*k];
      auto val1imag = evenVec[2*k+1];
      auto val2real = std::cos(4*M_PI*k/numVals)*oddVec[2*k] - std::sin(4*M_PI*k/numVals)*oddVec[2*k+1];
      auto val2imag = std::sin(4*M_PI*k/numVals)*oddVec[2*k] + std::cos(4*M_PI*k/numVals)*oddVec[2*k+1];
      inVec[2*k] = val1real + val2real;
      inVec[2*k+1] = val1imag + val2imag;
      inVec[2*k+(numVals/2)] = val1real - val2real;
      inVec[2*k+(numVals/2)+1]=  val1imag - val2imag;
    }
  }
}

// complex data type implementation.  Easier to debug.
void calcCTForwardFFTComplex( std::vector<std::complex<double> > & inVec)
{
  auto numVals = inVec.size();
  if( numVals == 1)
  {
    // trivial case of size 1
  }
  else
  {
    // split the input vector 
    std::vector<std::complex<double> > evenVec(numVals/2);
    std::vector<std::complex<double> > oddVec(numVals/2);
    for( auto j=0; j< (numVals/2); j++)
    {
      evenVec[j] = inVec[2*j];
      oddVec[j] = inVec[2*j+1];
    }
    calcCTForwardFFTComplex( evenVec );
    calcCTForwardFFTComplex( oddVec);
    // now combine even and odd elements.
    for( auto k = 0; k<(numVals/2); k++)
    {
      const std::complex<double> kernel(0.0,-2*M_PI*k/numVals);
      auto val1 = evenVec[k];
      auto val2 = std::exp( kernel )*oddVec[k];
      inVec[k] = val1 + val2;
      inVec[k+(numVals/2)] = val1 - val2;
    }
  }
}



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
// Purpose        : Unit tests for FFT interface class
//
// Special Notes  : This tests the FFT interface with more complicated 
//                  inputs than what was tested in the old unit-test code.
//
// Creator        : Richard Schiek, SNL, Parallel Computational Sciences
//
// Creation Date  : 2024/2022
//
//
//
//
//-------------------------------------------------------------------------

#include <gtest/gtest.h>
#include "Xyce_config.h"
#include <N_UTL_FFTInterface.hpp>
#include <math.h>
#include <iostream>
#include <sstream>

void ReadInTestData( const std::string aFileName,  const size_t expectedNumPoints,
  std::vector<double> & timeVec,
  std::vector<double> & functionVec,
  std::vector<double> & fftRealVec, 
  std::vector<double> & fftImaglVec,
  std::vector<double> & ifftRealVec )
{
  // streams for doing the reading and parsing
  std::ifstream testdatastream;
  std::string aLine;
  
  testdatastream.open(aFileName);
  std::getline( testdatastream, aLine);
  std::istringstream streamConverter1(aLine);
  double aNumber=0.0;
  while( streamConverter1 >> aNumber)
  {
    timeVec.push_back(aNumber);
  }
  EXPECT_EQ( timeVec.size(), expectedNumPoints);
  
  std::getline( testdatastream, aLine);
  std::istringstream streamConverter2(aLine);
  while( streamConverter2 >> aNumber)
  {
    functionVec.push_back(aNumber);
  }
  EXPECT_EQ( functionVec.size(), expectedNumPoints);
  std::getline( testdatastream, aLine);
  std::istringstream streamConverter3(aLine);
  while( streamConverter3 >> aNumber)
  {
    fftRealVec.push_back(aNumber);
  }
  EXPECT_EQ( fftRealVec.size(), expectedNumPoints);
  std::getline( testdatastream, aLine);
  std::istringstream streamConverter4(aLine);
  while( streamConverter4 >> aNumber)
  {
    fftImaglVec.push_back(aNumber);
  }
  EXPECT_EQ( fftImaglVec.size(), expectedNumPoints);
  std::getline( testdatastream, aLine);
  std::istringstream streamConverter5(aLine);
  while( streamConverter5 >> aNumber)
  {
    ifftRealVec.push_back(aNumber);
  }
  EXPECT_EQ( ifftRealVec.size(), expectedNumPoints);
  testdatastream.close();
}

TEST ( Power2, create)
{
  // use a power of 2 number of points 
  const int numPts = (int)(std::pow(2, 8));
  // create the interface  
  N_UTL_FFTInterface<std::vector<double> > * fftInterfacePtr = NULL;
  fftInterfacePtr = new  N_UTL_FFTInterface<std::vector<double> >( numPts );

  EXPECT_TRUE( fftInterfacePtr != NULL );
  delete fftInterfacePtr;
}

TEST ( NonPower2, create)
{
  // use a power of 2 number of points 
  const int numPts = (int)(std::pow(2, 8)) + 1;
  // create the interface  
  N_UTL_FFTInterface<std::vector<double> > * fftInterfacePtr = NULL;
  fftInterfacePtr = new  N_UTL_FFTInterface<std::vector<double> >( numPts );

  EXPECT_TRUE( fftInterfacePtr != NULL );
  delete fftInterfacePtr;
}


TEST ( Sine1FreqEven, FFT_IFFT)
{
  // data file name and expected size of input data
  const std::string filename("Sin1fEven.txt");
  const int expectedNumPoints = (int)(std::pow(2, 8));
  
  // get data for this test from file.
  // expected number of points is checked in the read function.
  std::vector<double> time, fxn, fxnFFTreal, fxnFFTimag, fxnIFFTreal;
  ReadInTestData(filename, expectedNumPoints, time, fxn, fxnFFTreal, fxnFFTimag, fxnIFFTreal);
  
  // create the interface  
  N_UTL_FFTInterface<std::vector<double> > * fftInterfacePtr = NULL;
  fftInterfacePtr = new  N_UTL_FFTInterface<std::vector<double> >( expectedNumPoints );
  EXPECT_TRUE( fftInterfacePtr != NULL );
  
  std::vector<double> outputSignal(expectedNumPoints+2, 0.0);
  std::vector<double> backSignal(expectedNumPoints, 0.0);
  fftInterfacePtr->calculateFFT( fxn, &outputSignal );
  // check the forward transform against input data
  for( auto i=0; i<expectedNumPoints; i=i+2)
  {
    //std::cerr << i << ": real=" << fxnFFTreal[i/2] << " == " << ": out " <<  outputSignal[i] << " AND imag=" << fxnFFTimag[1+i/2] << ", " << outputSignal[i+3] << std::endl;
    EXPECT_NEAR(fxnFFTreal[i/2], outputSignal[i], 1.0e-8 );
    EXPECT_NEAR(fxnFFTimag[1+i/2], outputSignal[i+3], 1.0e-8 );
  }
  fftInterfacePtr->calculateIFT( outputSignal, &backSignal );
  // check back transform against input data
  for( auto i=0; i<expectedNumPoints; i++)
  {
    EXPECT_NEAR(fxn[i], backSignal[i], 1.0e-8 );
  }
  
  delete fftInterfacePtr;
}


TEST ( Sine1FreqOdd, FFT_IFFT)
{
  // data file name and expected size of input data
  const std::string filename("Sin1fOdd.txt");
  const int expectedNumPoints = (int)(std::pow(2, 8))+17;
  
  // get data for this test from file.
  // expected number of points is checked in the read function.
  std::vector<double> time, fxn, fxnFFTreal, fxnFFTimag, fxnIFFTreal;
  ReadInTestData(filename, expectedNumPoints, time, fxn, fxnFFTreal, fxnFFTimag, fxnIFFTreal);
  
  // create the interface  
  N_UTL_FFTInterface<std::vector<double> > * fftInterfacePtr = NULL;
  fftInterfacePtr = new  N_UTL_FFTInterface<std::vector<double> >( expectedNumPoints );
  EXPECT_TRUE( fftInterfacePtr != NULL );
  
  std::vector<double> outputSignal(expectedNumPoints+1, 0.0);
  std::vector<double> backSignal(expectedNumPoints, 0.0);
  fftInterfacePtr->calculateFFT( fxn, &outputSignal );
  // check the forward transform against input data
  for( auto i=0; i<(expectedNumPoints-1); i=i+2)
  {
    //std::cerr << i << ": real=" << fxnFFTreal[i/2] << " == " << ": out " <<  outputSignal[i] << " AND imag=" << fxnFFTimag[1+i/2] << ", " << outputSignal[i+3] << std::endl;
    EXPECT_NEAR(fxnFFTreal[i/2], outputSignal[i], 5.0e-8 );
    EXPECT_NEAR(fxnFFTimag[1+i/2], outputSignal[i+3], 5.0e-8 );
  }
  fftInterfacePtr->calculateIFT( outputSignal, &backSignal );
  // check back transform against input data
  for( auto i=0; i<expectedNumPoints; i++)
  {
    EXPECT_NEAR(fxn[i], backSignal[i], 1.0e-8 );
  }
  
  delete fftInterfacePtr;
}



TEST ( Sine2FreqEven, FFT_IFFT)
{
  // data file name and expected size of input data
  const std::string filename("Sin2fEven.txt");
  const int expectedNumPoints = (int)(std::pow(2, 8));
  
  // get data for this test from file.
  // expected number of points is checked in the read function.
  std::vector<double> time, fxn, fxnFFTreal, fxnFFTimag, fxnIFFTreal;
  ReadInTestData(filename, expectedNumPoints, time, fxn, fxnFFTreal, fxnFFTimag, fxnIFFTreal);
  
  // create the interface  
  N_UTL_FFTInterface<std::vector<double> > * fftInterfacePtr = NULL;
  fftInterfacePtr = new  N_UTL_FFTInterface<std::vector<double> >( expectedNumPoints );
  EXPECT_TRUE( fftInterfacePtr != NULL );
  
  std::vector<double> outputSignal(expectedNumPoints+2, 0.0);
  std::vector<double> backSignal(expectedNumPoints, 0.0);
  fftInterfacePtr->calculateFFT( fxn, &outputSignal );
  // check the forward transform against input data
  for( auto i=0; i<expectedNumPoints; i=i+2)
  {
    //std::cerr << i << ": real=" << fxnFFTreal[i/2] << " == " << ": out " <<  outputSignal[i] << " AND imag=" << fxnFFTimag[1+i/2] << ", " << outputSignal[i+3] << std::endl;
    EXPECT_NEAR(fxnFFTreal[i/2], outputSignal[i], 1.0e-8 );
    EXPECT_NEAR(fxnFFTimag[1+i/2], outputSignal[i+3], 1.0e-8 );
  }
  fftInterfacePtr->calculateIFT( outputSignal, &backSignal );
  // check back transform against input data
  for( auto i=0; i<expectedNumPoints; i++)
  {
    EXPECT_NEAR(fxn[i], backSignal[i], 1.0e-8 );
  }
  
  delete fftInterfacePtr;
}

TEST ( Sine2FreqOdd, FFT_IFFT)
{
  // data file name and expected size of input data
  const std::string filename("Sin2fOdd.txt");
  const int expectedNumPoints = (int)(std::pow(2, 8))+17;
  
  // get data for this test from file.
  // expected number of points is checked in the read function.
  std::vector<double> time, fxn, fxnFFTreal, fxnFFTimag, fxnIFFTreal;
  ReadInTestData(filename, expectedNumPoints, time, fxn, fxnFFTreal, fxnFFTimag, fxnIFFTreal);
  
  // create the interface  
  N_UTL_FFTInterface<std::vector<double> > * fftInterfacePtr = NULL;
  fftInterfacePtr = new  N_UTL_FFTInterface<std::vector<double> >( expectedNumPoints );
  EXPECT_TRUE( fftInterfacePtr != NULL );
  
  std::vector<double> outputSignal(expectedNumPoints+1, 0.0);
  std::vector<double> backSignal(expectedNumPoints, 0.0);
  fftInterfacePtr->calculateFFT( fxn, &outputSignal );
  // check the forward transform against input data
  for( auto i=0; i<(expectedNumPoints-1); i=i+2)
  {
    //std::cerr << i << ": real=" << fxnFFTreal[i/2] << " == " << ": out " <<  outputSignal[i] << " AND imag=" << fxnFFTimag[1+i/2] << ", " << outputSignal[i+3] << std::endl;
    EXPECT_NEAR(fxnFFTreal[i/2], outputSignal[i], 5.0e-8 );
    EXPECT_NEAR(fxnFFTimag[1+i/2], outputSignal[i+3], 5.0e-8 );
  }
  fftInterfacePtr->calculateIFT( outputSignal, &backSignal );
  // check back transform against input data
  for( auto i=0; i<expectedNumPoints; i++)
  {
    EXPECT_NEAR(fxn[i], backSignal[i], 1.0e-8 );
  }
  
  delete fftInterfacePtr;
}

TEST ( Sine3FreqEven, FFT_IFFT)
{
  // data file name and expected size of input data
  const std::string filename("Sin3fEven.txt");
  const int expectedNumPoints = (int)(std::pow(2, 8));
  
  // get data for this test from file.
  // expected number of points is checked in the read function.
  std::vector<double> time, fxn, fxnFFTreal, fxnFFTimag, fxnIFFTreal;
  ReadInTestData(filename, expectedNumPoints, time, fxn, fxnFFTreal, fxnFFTimag, fxnIFFTreal);
  
  // create the interface  
  N_UTL_FFTInterface<std::vector<double> > * fftInterfacePtr = NULL;
  fftInterfacePtr = new  N_UTL_FFTInterface<std::vector<double> >( expectedNumPoints );
  EXPECT_TRUE( fftInterfacePtr != NULL );
  
  std::vector<double> outputSignal(expectedNumPoints+2, 0.0);
  std::vector<double> backSignal(expectedNumPoints, 0.0);
  fftInterfacePtr->calculateFFT( fxn, &outputSignal );
  // check the forward transform against input data
  for( auto i=0; i<expectedNumPoints; i=i+2)
  {
    //std::cerr << i << ": real=" << fxnFFTreal[i/2] << " == " << ": out " <<  outputSignal[i] << " AND imag=" << fxnFFTimag[1+i/2] << ", " << outputSignal[i+3] << std::endl;
    EXPECT_NEAR(fxnFFTreal[i/2], outputSignal[i], 1.0e-8 );
    EXPECT_NEAR(fxnFFTimag[1+i/2], outputSignal[i+3], 1.0e-8 );
  }
  fftInterfacePtr->calculateIFT( outputSignal, &backSignal );
  // check back transform against input data
  for( auto i=0; i<expectedNumPoints; i++)
  {
    EXPECT_NEAR(fxn[i], backSignal[i], 1.0e-8 );
  }
  
  delete fftInterfacePtr;
}

TEST ( Sine3FreqOdd, FFT_IFFT)
{
  // data file name and expected size of input data
  const std::string filename("Sin3fOdd.txt");
  const int expectedNumPoints = (int)(std::pow(2, 8))+17;
  
  // get data for this test from file.
  // expected number of points is checked in the read function.
  std::vector<double> time, fxn, fxnFFTreal, fxnFFTimag, fxnIFFTreal;
  ReadInTestData(filename, expectedNumPoints, time, fxn, fxnFFTreal, fxnFFTimag, fxnIFFTreal);
  
  // create the interface  
  N_UTL_FFTInterface<std::vector<double> > * fftInterfacePtr = NULL;
  fftInterfacePtr = new  N_UTL_FFTInterface<std::vector<double> >( expectedNumPoints );
  EXPECT_TRUE( fftInterfacePtr != NULL );
  
  std::vector<double> outputSignal(expectedNumPoints+1, 0.0);
  std::vector<double> backSignal(expectedNumPoints, 0.0);
  fftInterfacePtr->calculateFFT( fxn, &outputSignal );
  // check the forward transform against input data
  for( auto i=0; i<(expectedNumPoints-1); i=i+2)
  {
    //std::cerr << i << ": real=" << fxnFFTreal[i/2] << " == " << ": out " <<  outputSignal[i] << " AND imag=" << fxnFFTimag[1+i/2] << ", " << outputSignal[i+3] << std::endl;
    EXPECT_NEAR(fxnFFTreal[i/2], outputSignal[i], 5.0e-8 );
    EXPECT_NEAR(fxnFFTimag[1+i/2], outputSignal[i+3], 5.0e-8 );
  }
  fftInterfacePtr->calculateIFT( outputSignal, &backSignal );
  // check back transform against input data
  for( auto i=0; i<expectedNumPoints; i++)
  {
    EXPECT_NEAR(fxn[i], backSignal[i], 1.0e-8 );
  }
  
  delete fftInterfacePtr;
}


TEST ( Gauss1Even, FFT_IFFT)
{
  // data file name and expected size of input data
  const std::string filename("Gauss1Even.txt");
  const int expectedNumPoints = (int)(std::pow(2, 8));
  
  // get data for this test from file.
  // expected number of points is checked in the read function.
  std::vector<double> time, fxn, fxnFFTreal, fxnFFTimag, fxnIFFTreal;
  ReadInTestData(filename, expectedNumPoints, time, fxn, fxnFFTreal, fxnFFTimag, fxnIFFTreal);
  
  // create the interface  
  N_UTL_FFTInterface<std::vector<double> > * fftInterfacePtr = NULL;
  fftInterfacePtr = new  N_UTL_FFTInterface<std::vector<double> >( expectedNumPoints );
  EXPECT_TRUE( fftInterfacePtr != NULL );
  
  std::vector<double> outputSignal(expectedNumPoints+2, 0.0);
  std::vector<double> backSignal(expectedNumPoints, 0.0);
  fftInterfacePtr->calculateFFT( fxn, &outputSignal );
  // check the forward transform against input data
  for( auto i=0; i<expectedNumPoints; i=i+2)
  {
    //std::cerr << i << ": real=" << fxnFFTreal[i/2] << " == " << ": out " <<  outputSignal[i] << " AND imag=" << fxnFFTimag[1+i/2] << ", " << outputSignal[i+3] << std::endl;
    EXPECT_NEAR(fxnFFTreal[i/2], outputSignal[i], 1.0e-8 );
    EXPECT_NEAR(fxnFFTimag[1+i/2], outputSignal[i+3], 1.0e-8 );
  }
  fftInterfacePtr->calculateIFT( outputSignal, &backSignal );
  // check back transform against input data
  for( auto i=0; i<expectedNumPoints; i++)
  {
    EXPECT_NEAR(fxn[i], backSignal[i], 1.0e-8 );
  }
  
  delete fftInterfacePtr;
}

TEST ( Gauss1Odd, FFT_IFFT)
{
  // data file name and expected size of input data
  const std::string filename("Gauss1Odd.txt");
  const int expectedNumPoints = (int)(std::pow(2, 8))+17;
  
  // get data for this test from file.
  // expected number of points is checked in the read function.
  std::vector<double> time, fxn, fxnFFTreal, fxnFFTimag, fxnIFFTreal;
  ReadInTestData(filename, expectedNumPoints, time, fxn, fxnFFTreal, fxnFFTimag, fxnIFFTreal);
  
  // create the interface  
  N_UTL_FFTInterface<std::vector<double> > * fftInterfacePtr = NULL;
  fftInterfacePtr = new  N_UTL_FFTInterface<std::vector<double> >( expectedNumPoints );
  EXPECT_TRUE( fftInterfacePtr != NULL );
  
  std::vector<double> outputSignal(expectedNumPoints+1, 0.0);
  std::vector<double> backSignal(expectedNumPoints, 0.0);
  fftInterfacePtr->calculateFFT( fxn, &outputSignal );
  // check the forward transform against input data
  for( auto i=0; i<(expectedNumPoints-1); i=i+2)
  {
    //std::cerr << i << ": real=" << fxnFFTreal[i/2] << " == " << ": out " <<  outputSignal[i] << " AND imag=" << fxnFFTimag[1+i/2] << ", " << outputSignal[i+3] << std::endl;
    EXPECT_NEAR(fxnFFTreal[i/2], outputSignal[i], 1.0e-7 );
    EXPECT_NEAR(fxnFFTimag[1+i/2], outputSignal[i+3], 1.0e-7 );
  }
  fftInterfacePtr->calculateIFT( outputSignal, &backSignal );
  // check back transform against input data
  for( auto i=0; i<expectedNumPoints; i++)
  {
    EXPECT_NEAR(fxn[i], backSignal[i], 1.0e-8 );
  }
  
  delete fftInterfacePtr;
}


TEST ( Gauss2Even, FFT_IFFT)
{
  // data file name and expected size of input data
  const std::string filename("Gauss2Even.txt");
  const int expectedNumPoints = (int)(std::pow(2, 8));
  
  // get data for this test from file.
  // expected number of points is checked in the read function.
  std::vector<double> time, fxn, fxnFFTreal, fxnFFTimag, fxnIFFTreal;
  ReadInTestData(filename, expectedNumPoints, time, fxn, fxnFFTreal, fxnFFTimag, fxnIFFTreal);
  
  // create the interface  
  N_UTL_FFTInterface<std::vector<double> > * fftInterfacePtr = NULL;
  fftInterfacePtr = new  N_UTL_FFTInterface<std::vector<double> >( expectedNumPoints );
  EXPECT_TRUE( fftInterfacePtr != NULL );
  
  std::vector<double> outputSignal(expectedNumPoints+2, 0.0);
  std::vector<double> backSignal(expectedNumPoints, 0.0);
  fftInterfacePtr->calculateFFT( fxn, &outputSignal );
  // check the forward transform against input data
  for( auto i=0; i<expectedNumPoints; i=i+2)
  {
    //std::cerr << i << ": real=" << fxnFFTreal[i/2] << " == " << ": out " <<  outputSignal[i] << " AND imag=" << fxnFFTimag[1+i/2] << ", " << outputSignal[i+3] << std::endl;
    EXPECT_NEAR(fxnFFTreal[i/2], outputSignal[i], 1.0e-7 );
    EXPECT_NEAR(fxnFFTimag[1+i/2], outputSignal[i+3], 1.0e-7 );
  }
  fftInterfacePtr->calculateIFT( outputSignal, &backSignal );
  // check back transform against input data
  for( auto i=0; i<expectedNumPoints; i++)
  {
    EXPECT_NEAR(fxn[i], backSignal[i], 1.0e-8 );
  }
  
  delete fftInterfacePtr;
}

TEST ( Gauss2Odd, FFT_IFFT)
{
  // data file name and expected size of input data
  const std::string filename("Gauss2Odd.txt");
  const int expectedNumPoints = (int)(std::pow(2, 8))+17;
  
  // get data for this test from file.
  // expected number of points is checked in the read function.
  std::vector<double> time, fxn, fxnFFTreal, fxnFFTimag, fxnIFFTreal;
  ReadInTestData(filename, expectedNumPoints, time, fxn, fxnFFTreal, fxnFFTimag, fxnIFFTreal);
  
  // create the interface  
  N_UTL_FFTInterface<std::vector<double> > * fftInterfacePtr = NULL;
  fftInterfacePtr = new  N_UTL_FFTInterface<std::vector<double> >( expectedNumPoints );
  EXPECT_TRUE( fftInterfacePtr != NULL );
  
  std::vector<double> outputSignal(expectedNumPoints+1, 0.0);
  std::vector<double> backSignal(expectedNumPoints, 0.0);
  fftInterfacePtr->calculateFFT( fxn, &outputSignal );
  // check the forward transform against input data
  for( auto i=0; i<(expectedNumPoints-1); i=i+2)
  {
    //std::cerr << i << ": real=" << fxnFFTreal[i/2] << " == " << ": out " <<  outputSignal[i] << " AND imag=" << fxnFFTimag[1+i/2] << ", " << outputSignal[i+3] << std::endl;
    EXPECT_NEAR(fxnFFTreal[i/2], outputSignal[i], 1.0e-7 );
    EXPECT_NEAR(fxnFFTimag[1+i/2], outputSignal[i+3], 1.0e-7 );
  }
  fftInterfacePtr->calculateIFT( outputSignal, &backSignal );
  // check back transform against input data
  for( auto i=0; i<expectedNumPoints; i++)
  {
    EXPECT_NEAR(fxn[i], backSignal[i], 1.0e-8 );
  }
  
  delete fftInterfacePtr;
}

TEST ( Step1Even, FFT_IFFT)
{
  // data file name and expected size of input data
  const std::string filename("Step1Even.txt");
  const int expectedNumPoints = (int)(std::pow(2, 8));
  
  // get data for this test from file.
  // expected number of points is checked in the read function.
  std::vector<double> time, fxn, fxnFFTreal, fxnFFTimag, fxnIFFTreal;
  ReadInTestData(filename, expectedNumPoints, time, fxn, fxnFFTreal, fxnFFTimag, fxnIFFTreal);
  
  // create the interface  
  N_UTL_FFTInterface<std::vector<double> > * fftInterfacePtr = NULL;
  fftInterfacePtr = new  N_UTL_FFTInterface<std::vector<double> >( expectedNumPoints );
  EXPECT_TRUE( fftInterfacePtr != NULL );
  
  std::vector<double> outputSignal(expectedNumPoints+2, 0.0);
  std::vector<double> backSignal(expectedNumPoints, 0.0);
  fftInterfacePtr->calculateFFT( fxn, &outputSignal );
  // check the forward transform against input data
  for( auto i=0; i<expectedNumPoints; i=i+2)
  {
    //std::cerr << i << ": real=" << fxnFFTreal[i/2] << " == " << ": out " <<  outputSignal[i] << " AND imag=" << fxnFFTimag[1+i/2] << ", " << outputSignal[i+3] << std::endl;
    EXPECT_NEAR(fxnFFTreal[i/2], outputSignal[i], 1.0e-8 );
    EXPECT_NEAR(fxnFFTimag[1+i/2], outputSignal[i+3], 1.0e-8 );
  }
  fftInterfacePtr->calculateIFT( outputSignal, &backSignal );
  // check back transform against input data
  for( auto i=0; i<expectedNumPoints; i++)
  {
    EXPECT_NEAR(fxn[i], backSignal[i], 1.0e-8 );
  }
  
  delete fftInterfacePtr;
}

TEST ( Step1Odd, FFT_IFFT)
{
  // data file name and expected size of input data
  const std::string filename("Step1Odd.txt");
  const int expectedNumPoints = (int)(std::pow(2, 8))+17;
  
  // get data for this test from file.
  // expected number of points is checked in the read function.
  std::vector<double> time, fxn, fxnFFTreal, fxnFFTimag, fxnIFFTreal;
  ReadInTestData(filename, expectedNumPoints, time, fxn, fxnFFTreal, fxnFFTimag, fxnIFFTreal);
  
  // create the interface  
  N_UTL_FFTInterface<std::vector<double> > * fftInterfacePtr = NULL;
  fftInterfacePtr = new  N_UTL_FFTInterface<std::vector<double> >( expectedNumPoints );
  EXPECT_TRUE( fftInterfacePtr != NULL );
  
  std::vector<double> outputSignal(expectedNumPoints+1, 0.0);
  std::vector<double> backSignal(expectedNumPoints, 0.0);
  fftInterfacePtr->calculateFFT( fxn, &outputSignal );
  // check the forward transform against input data
  for( auto i=0; i<(expectedNumPoints-1); i=i+2)
  {
    //std::cerr << i << ": real=" << fxnFFTreal[i/2] << " == " << ": out " <<  outputSignal[i] << " AND imag=" << fxnFFTimag[1+i/2] << ", " << outputSignal[i+3] << std::endl;
    EXPECT_NEAR(fxnFFTreal[i/2], outputSignal[i], 1.0e-7 );
    EXPECT_NEAR(fxnFFTimag[1+i/2], outputSignal[i+3], 1.0e-7 );
  }
  fftInterfacePtr->calculateIFT( outputSignal, &backSignal );
  // check back transform against input data
  for( auto i=0; i<expectedNumPoints; i++)
  {
    EXPECT_NEAR(fxn[i], backSignal[i], 1.0e-8 );
  }
  
  delete fftInterfacePtr;
}


//-------------------------------------------------------------------------------
int main (int argc, char **argv)
{
  testing::InitGoogleTest(&argc, argv);

  return RUN_ALL_TESTS();
}
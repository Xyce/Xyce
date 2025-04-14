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

//
// test the N_UTL_FFTInterface class
//

#include <N_UTL_FFTInterface.hpp>
#include <N_UTL_fwd.h>
#include <N_PDS_ParMap.h>
#include <N_PDS_Comm.h>
#include <N_PDS_ParHelpers.h>
#include <N_LAS_EpetraVector.h>
#include <N_LAS_BlockVector.h>

#include <Epetra_SerialComm.h>
#include <Epetra_Map.h>
#include <Epetra_BlockMap.h>

#include <iostream>
#include <vector>
#include <N_UTL_Math.h>

int main(int argc, char* argv[])
{
  //
  // This first part of the code tests out making a signal x(t) and then
  // taking its fourier transform fft( x(t) ) and then its back 
  // transform ift( fft( x(t) ) ).  This involves running the FFT interface
  // with a single input vector (either a vector<double> or an N_LAS_Vector
  // object.  In the second half of this code we try doing multiple FFTs at
  // the same time using an N_LAS_BlockVector object
  //
  
  // create arrays to hold some data for testing
  //int numPts = 128;
  int numPts = 11;
  int lengthTransformedSignal = numPts + 2;
  if( numPts % 2 )
  {
    // numPts was odd so the actual number of elements in the transformed signal is numPts+1
    lengthTransformedSignal = numPts + 1;
  }
  std::cout << "numPts = " << numPts << std::endl;
  std::cout << "lengthTransformedSignal = " << lengthTransformedSignal << std::endl;
  
  double timeStart = 0.0;
  double timeStop = 1;
  double deltaTime = (timeStop - timeStart)/numPts;
  double freqDelta = 1.0 / (timeStop - timeStart);
  
  // these are to test the vector of double accecess to FFT interface 
  std::vector<double> time(numPts, 0.0);
  std::vector<double> inputSignal(numPts, 0.0);
  std::vector<double> outputSignal(lengthTransformedSignal, 0.0);
  std::vector<double> backSignal(numPts, 0.0);
  
  // these is to make Xyce::Linear::MultiVector objects to test that part of
  // the interface since Xyce::Linear::Vector derives from MultiVector, one
  // can use this same interface for Xyce::Linear::Vector objects
  Teuchos::RCP<Xyce::Parallel::Communicator> pdsComm = Teuchos::rcp( Xyce::Parallel::createPDSComm( argc, argv ) );
  std::vector<int> realLbMap(numPts, 0.0);
  std::vector<int> complexLbMap(lengthTransformedSignal, 0.0);
  Teuchos::RCP<Xyce::Parallel::ParMap> parMapForReal = 
    Teuchos::rcp( Xyce::Parallel::createPDSParMap( numPts, numPts, realLbMap, 0, *pdsComm.get() ) );
  Teuchos::RCP<Xyce::Parallel::ParMap> parMapForComplex = 
    Teuchos::rcp( Xyce::Parallel::createPDSParMap( lengthTransformedSignal, lengthTransformedSignal, complexLbMap, 0, *pdsComm.get() ) );
  
  Xyce::Linear::EpetraVector timeMV( *parMapForReal );
  Xyce::Linear::EpetraVector inputSignalMV( *parMapForReal );
  Xyce::Linear::EpetraVector outputSignalMV( *parMapForComplex );
  Xyce::Linear::EpetraVector backSignalMV( *parMapForReal );
  
  // Fill in the vectors
  for(int i=0; i<numPts; i++)
  {
    time[i] = timeStart + i * deltaTime;
    
    inputSignal[i] = std::sin( 2.0 * M_PI * 1 * time[i] ) + time[i]*time[i];      }

  // create the interface  
  N_UTL_FFTInterface<std::vector<double> > myTransform( numPts );
  
  // try taking an FFT
  myTransform.calculateFFT( inputSignal, &outputSignal );
  myTransform.calculateIFT( outputSignal, &backSignal );
  
  std::cout.precision(6);
  Xyce::dout() << "time\tInputSignal\tBackSignal" << std::endl;
  for(int i=0; i<numPts; i++)
  {
    Xyce::dout() << time[i] << "\t" << inputSignal[i] << "\t" << backSignal[i] << std::endl;
     
  }
  Xyce::dout() << "Frequency\tReal + Imaginary" << std::endl;
  for (int i=0 ; i<lengthTransformedSignal/2  ; ++i) {
      Xyce::dout() << i*freqDelta << "\t" << outputSignal[2*i] << " + "
        << outputSignal[2*i+1] << "i" << std::endl;
  }
  
  return 0;
}

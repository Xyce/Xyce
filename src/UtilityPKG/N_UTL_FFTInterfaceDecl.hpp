//-------------------------------------------------------------------------
//   Copyright 2002-2019 National Technology & Engineering Solutions of
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
// Purpose        : This class acts as an abstract interface to an FFT library
//                  for FFT and IFT calculations.  This class should isolate
//                  Xyce from the specifics of a given FFT library so 
//                  that multiple libraries can be used.  
//
// Special Notes  : 
//
// Creator        : Heidi Thornquist
//
// Creation Date  : 11/11/10
//
//-------------------------------------------------------------------------
#ifndef N_UTL_FFTINTERFACE_DECL_HPP
#define N_UTL_FFTINTERFACE_DECL_HPP


// ---------- Standard Includes ----------

#include <N_UTL_DFTInterfaceDecl.hpp>

// ----------   Other Includes   ----------

#include <Teuchos_RCP.hpp>

// ---------- Structure definitions ----------

//-----------------------------------------------------------------------------
// Class         : FFTInterfaceDecl
// Purpose       : This class acts as an abstract interface to an FFT library
//                 for FFT and IFT calculations.  This class should isolate
//                 Xyce from the specifics of a given FFT library so 
//                 that multiple libraries can be used.  
// Special Notes :
// Creator       : Heidi Thornquist
// Creation Date : 11/11/10
// Last Modified : 11/11/10
//-----------------------------------------------------------------------------
template<typename VectorType>
class N_UTL_FFTInterfaceDecl : public N_UTL_DFTInterfaceDecl<VectorType>
{
  public:
    // Basic constructor without passing in vector structures, which need to be registered later or
    // passed into the calculate[FFT/IFT] methods.
    N_UTL_FFTInterfaceDecl( int length, int numSignals=1, int reqStride=0, bool overwrite=false )
      { signalLength_ = length; numberSignals_ = numSignals; stride_ = reqStride; overwrite_ = overwrite; }
    
    // Basic destructor 
    virtual ~N_UTL_FFTInterfaceDecl() {};

    // Return signal length.
    virtual int getSignalLength()
    { return signalLength_; }

    virtual int getScalar()
    { return signalLength_; }

    virtual void setSignalLength(int  length)
    { signalLength_ = length;}

  protected:
    // this is the length of the real array.  The complex result will be signalLength+2 long
    // if signalLength is even and signalLength+1 if signalLength is odd
    int signalLength_;

    // this is the number of signals on which we will take an fft/ift
    int numberSignals_;

    // If the signals are grouped by blocks at the same time (say x0, x1 ... xn at t0) and
    // then (x0, x1 ... xn at t1). Then stride is the spacing from one x0 at t0 to the next
    // x0 at t1.  This lets one take ffts/ifts of data that is blocked by time
    int stride_;

    // Whether the input and output vectors should be expected to be the same, save space if possible.
    bool overwrite_;
};

#endif


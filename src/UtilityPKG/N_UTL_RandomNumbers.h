//-------------------------------------------------------------------------
//   Copyright 2002-2022 National Technology & Engineering Solutions of
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
// Purpose        : Provide a class with methods to generate random numbers
//                  with specific properties
//
// Special Notes  :
//
// Creator        : Tom Russo, SNL
//
// Creation Date  :  5/28/2014
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef N_UTL_RandomNumbers_H
#define N_UTL_RandomNumbers_H

namespace Xyce {
namespace Util {

class RandomNumbers
{
 public:
  RandomNumbers(long seed=0, bool output=true)   
    : randInitialized_(false),
      useLastGaussian_(false),
      //poissonInitialized_(false),
      xPoisson_(-1.0),
      pPoisson_(-1.0),
      sPoisson_(-1.0)
    {
      seedRandom(seed,output);
    }

  double uniformRandom();
  double gaussianRandom(double mu, double sigma);
  int poissonRandom( double lambda );
  void seedRandom(long seed, bool output=true);
  int uniformRandomInt(int min=0, int max=1);

 private:
  bool randInitialized_;
#ifndef HAVE_DRAND48
  // These bits of data only needed if we have no good random number generator
  double maxran_;
  int randBuf_[98];
  double randY_;
#endif

  // Internals for use by gaussian random number generator
  double ySaveGaussian_;
  bool useLastGaussian_;
  // these are used in generating a Poisson random number from the 
  // uniform random number generator
  //bool poissonInitialized_;
  double xPoisson_;
  double pPoisson_;
  double sPoisson_;
};
}
}
#endif //N_UTL_RandomNumbers_H

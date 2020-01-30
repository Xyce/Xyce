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

//-------------------------------------------------------------------------
//
// Purpose       :
//
// Special Notes :
//
// Creator       : Tom Russo, SNL
//
// Creation Date : 5/28/2014
//
//
//
//
//-------------------------------------------------------------------------

#include <Xyce_config.h>

#include <cstdlib>
#include <N_UTL_Math.h>
#include <time.h>

#if 1
#include<iostream>
#endif

#include <N_UTL_fwd.h>
#include <N_UTL_RandomNumbers.h>

namespace Xyce {
namespace Util {
//-----------------------------------------------------------------------------
// Function      : Xyce::Util::RandomNumbers::gaussianRandom
// Purpose       : Provide standardized, portable source of high-quality
//                 random numbers, normally distributed with mean mu and
//                 standard deviation sigma.
// Special Notes : Method is a variant of the Box-Muller transformation.
//                 A pair of random numbers from a uniform distribution is
//                 selected such that they can be the coordinates of a
//                 unit vector.  The magnitude and phase of this vector
//                 are used in the Box-Muller transformation to create a 
//                 set of values that are normally distributed instead of
//                 uniformly distributed.
// Creator       : Tom Russo
// Creation Date : 05/27/2014
//-----------------------------------------------------------------------------
double RandomNumbers::gaussianRandom(double mu, double sigma)
{
  double x1, x2, w, y1;
  
  if (useLastGaussian_)   // use value from previous call
  {
    y1 = ySaveGaussian_;
    useLastGaussian_ = false;
  }
  else
  {
    do 
    {
      x1 = 2.0 * uniformRandom() - 1.0;
      x2 = 2.0 * uniformRandom() - 1.0;
      w = x1 * x1 + x2 * x2;
    } while ( w >= 1.0 );
    
    w = sqrt( (-2.0 * log( w ) ) / w );
    y1 = x1 * w;
    ySaveGaussian_ = x2 * w;
    useLastGaussian_ = true;
  }
  return( mu + y1 * sigma );
}

//-----------------------------------------------------------------------------
// Function      : Xyce::Util::RandomNumbers::uniformRandom
// Purpose       : Provide standardized, portable source of high-quality
//                 random numbers, uniformly distributed on [0,1].
// Special Notes : When drand48 is available, we use it.  When only the
//                 very poor rand() function is available, it is used to
//                 create a pool of random numbers that are further reshuffled
//                 to improve the randomness of the sequence.
// Creator       : Tom Russo
// Creation Date : 05/21/2014
//-----------------------------------------------------------------------------
double RandomNumbers::uniformRandom ()
{

  double res;
  
#ifdef HAVE_DRAND48
  // The rand48 package is available on most Unix-like systems, and provides
  // a high-quality random numbers.  Just use it.
  res=drand48();
#else
  
  // Otherwise, all we can count on having is the generic "rand"
  // function, which is well known to be horrid.  Let's use it to create a 
  // better choice of random numbers.
  double dum;
  int j;

  if (!randInitialized_)
  {
    randInitialized_=true;
    maxran_=RAND_MAX+1.0;
    for (j=0;j<98;++j) dum=rand();
    for (j=0;j<98;++j) randBuf_[j]=rand();
    randY_=rand();
  }

  //  j will be a number 0-97
  j=98.0*randY_/maxran_;
  randY_=randBuf_[j];
  randBuf_[j]=rand();
  res = randY_/maxran_;
#endif

  return res;
}

//-----------------------------------------------------------------------------
// Function      : Xyce::Util::RandomNumbers::seedRandom
// Purpose       : Initialize a random number generator
// Special Notes : 
// Creator       : Tom Russo
// Creation Date : 05/21/2014
//-----------------------------------------------------------------------------
void RandomNumbers::seedRandom(long seed, bool output)
{
  if (seed==0)
  {
    seed=time(NULL);
  }

  if (output)
  {
    Xyce::lout() << "Seeding random number generator with " << seed << std::endl;
  }
  
#ifdef HAVE_DRAND48
  srand48(seed);
#else
  srand(seed);
#endif

}


//-----------------------------------------------------------------------------
// Function      : Xyce::Util::RandomNumbers::poissonRandom
// Purpose       : Return a Poisson distributed random number
// Special Notes : Based on: Devroye, Luc (1986). "Discrete Univariate Distributions". 
//                 Non-Uniform Random Variate Generation. New York: Springer-Verlag. p. 505.
// Creator       : Rich Schiek
// Creation Date : 04/10/2015
//-----------------------------------------------------------------------------
int RandomNumbers::poissonRandom( double lambda )
{
  //if( !poissonInitialized_ )
  //{
  //  poissonInitialized_=true;
  xPoisson_=0.0;
  pPoisson_ = exp( -lambda );
  sPoisson_ = pPoisson_;
  //}
  double u = uniformRandom();
  while (u > sPoisson_)
  {
    xPoisson_++;
    pPoisson_ *= lambda / xPoisson_;
    sPoisson_ += pPoisson_;  
  }
  return xPoisson_;
}

//-----------------------------------------------------------------------------
// Function      : Xyce::Util::RandomNumbers::uniformRandomInt
// Purpose       : Computes a uniformly distributed random integer 
//                 number between min and max.
//
//                 This function uses this formula:
//
//                     min + (rand() % static_cast<int>(max - min + 1));
//
// (from https://stackoverflow.com/questions/5008804/generating-random-integer-from-a-range )
//
// rand() not a very good random number generator, and to mitigate this, I am calling the
// "uniformRandom" function instead, as it uses a modified rand() when drand48() is not 
// available.
//
// This is not a great solution but will suffice until we switch to c++11 and start using <random>
//
// Special Notes : 
// Creator       : Eric Keiter
// Creation Date : 8/24/2018
//-----------------------------------------------------------------------------
int RandomNumbers::uniformRandomInt (int min, int max)
{
  int tmp = static_cast<int> (uniformRandom ()*RAND_MAX);
  return (min + (tmp % static_cast<int>(max - min + 1)));
}

} // Util
} // Xyce

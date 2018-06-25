/*
 * CheckForLinuxExceptions.cxx
 *
 *  Created on: Apr 17, 2018
 *      Author: asgibso
 */

#include <fenv.h>

using namespace std;

int main(int argc, char** argv)
{
  feenableexcept(FE_DIVBYZERO | FE_INVALID);
  return 0;
}



/*
 * CheckForIota.cpp
 *
 *  Created on: Apr 2, 2018
 *      Author: asgibso
 */

#include <vector>
#include <numeric>

int main(int argc, char** argv)
{
  std::vector<int> v(5, 0);
  std::iota(v.begin(), v.end(), 1);

  return 0;
}



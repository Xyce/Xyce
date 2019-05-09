#include <algorithm>
#include <numeric>
#include <vector>


int main ( void ) 
{ 
  std::vector<int> ipVec( 5, 0 );

  std::iota( ipVec.begin(), ipVec.end(), 1 );

  return 0; 
}

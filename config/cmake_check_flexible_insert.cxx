#include <vector>
#include <list>

int main ( void ) 
{ 
  std::list<int> ipList(1,0); 
  std::vector<int> ipVec; 

  ipVec.insert(ipVec.begin(),ipList.begin(),ipList.end()); 

  return 0; 
}

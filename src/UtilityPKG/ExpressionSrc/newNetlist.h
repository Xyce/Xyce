
#include<string>

#ifndef newNetlist_H
#define newNetlist_H

#include "netlistData.h"

namespace Xyce {
namespace Util {
void lexAndParseNetlist(std::string & netlistString, netlistData & nd);
}
}

#endif


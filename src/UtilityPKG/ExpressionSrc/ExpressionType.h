
#ifndef ExpressionType_H
#define ExpressionType_H

#include <algorithm>

typedef double basicType;

#ifdef USE_TYPE_DOUBLE
typedef double usedType ;
#else
typedef std::complex<basicType> usedType ;
#endif

#endif

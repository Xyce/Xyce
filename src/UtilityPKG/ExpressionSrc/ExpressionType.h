
#ifndef ExpressionType_H
#define ExpressionType_H

#include <algorithm>

typedef double basicType;

#ifdef USE_TYPE_COMPLEX
typedef std::complex<basicType> usedType ;
#else
typedef double usedType ;
#endif

#endif

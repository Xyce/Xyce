
#ifndef ExpressionType_H
#define ExpressionType_H

#include <algorithm>

typedef double basicType;

#ifdef USE_TYPE_COMPLEX
typedef std::complex<basicType> usedType ;
#else
typedef double usedType ;
#endif

//typedef Sacado::Fad::SFad<basicType,1> fadType;
//typedef Sacado::Fad::DFad<basicType> fadType;
//typedef fadType usedType ;

#endif

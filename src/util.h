#ifndef _UTIL_H_
#define _UTIL_H_

// fix Microsoft VC++ oddities
#ifdef min // MS defines min as a macro, making it impossible to use the identifier for functions etc.
#undef min
#endif
#ifdef max // MS defines max as a macro, making it impossible to use the identifier for functions etc.
#undef max
#endif
#ifndef M_PI // MS fails to define M_PI in the standard math header file
#define M_PI 3.14159265358979323846
#endif


#include "SSEVector3.h"
#include <xmmintrin.h>

template<class T> inline T sqr(const T &x) { return x*x; }

template<class T> inline T max(const T &a1, const T &a2, const T &a3)
{ 
	if(a1>a2) 
		return max(a1,a3); 
	else 
		return max(a2,a3);
}

template<class T> inline T max(const T &a1, const T &a2)
{ 
	if(a1>a2) 
		return a1; 
	else 
		return a2; 
}

#endif
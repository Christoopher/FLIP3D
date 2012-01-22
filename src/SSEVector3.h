#ifndef _SSE_VECTOR3_H_
#define _SSE_VECTOR3_H_

#include <iostream>
#include <cmath>
#include <xmmintrin.h>
#include <intrin.h>
#include "util.h"


namespace SSE
{

	template<class T>
	inline T sqr(const T &x) { return x*x; }

struct vec3f
{
	float __declspec(align(16)) v[3];

	vec3f() 
	{
		v[0] = v[1] = v[2] = 0; 
	}

	vec3f(float v_) { v[0] = v[1] = v[2] = v_; 	}

	vec3f(float x, float y, float z) { v[0] = x; v[1] = y; v[2] = z; }

	vec3f(vec3f & source) { v[0] = source.v[0]; v[1] = source.v[1]; v[2] = source.v[2]; }

	float & operator[](int index) {return v[index]; }

	const vec3f & operator[](int index ) const { return v[index]; }

	vec3f operator+=(const vec3f & rhs)
	{
		__m128 * pRhs =(__m128 *) &rhs;
		__m128 * pV =(__m128 *) v;

		*pV = _mm_add_ps(*pV,*pRhs);		
		return *this;
	}

	vec3f operator-=(const vec3f & rhs)
	{
		__m128 * pRhs =(__m128 *) &rhs;
		__m128 * pV =(__m128 *) v;

		*pV = _mm_sub_ps(*pV,*pRhs);		
		return *this;
	}

	vec3f operator*=(const float & scalar)
	{
		 __m128 mScalar = _mm_set_ps1(scalar);
		__m128 * pV =(__m128 *) v;

		*pV = _mm_mul_ps(*pV,mScalar);		
		return *this;
	}

	vec3f operator/=(const float & scalar)
	{
		__m128 mScalar = _mm_set_ps1(scalar);
		__m128 * pV =(__m128 *) v;

		*pV = _mm_div_ps(*pV,mScalar);		
		return *this;
	}
}; 


inline float mag2(const vec3f &a)
{ 
	
	__m128 ma = _mm_setr_ps(a.v[0],a.v[1],a.v[2],0.0f);
	float sum[4];
	__m128 * msum = (__m128 *) sum;
	*msum = _mm_mul_ps(ma, ma); // sum = { v0*a0, v1*a1, v2*a2, 0 }
	*msum = _mm_hadd_ps(*msum, *msum); // sum = { v0*a0 + v1*a1, v2*a2 + 0, v0*a0 + v1*a1, v2*a2 + 0}
	*msum = _mm_hadd_ps(*msum, *msum); // sum = { v0*a0 + v1*a1 + v2*a2 + 0 , ... }
	
	return sum[0];

	//return a.v[0]*a.v[0] + a.v[1]*a.v[1] + a.v[2]*a.v[2]; 
}

inline float mag(const vec3f &a)
{ return std::sqrtf(mag2(a)); }

inline float dist2(const vec3f &a, const vec3f &b)
{ return sqr(a.v[0]-b.v[0]) + sqr(a.v[1]-b.v[1]) + sqr(a.v[2]-b.v[2]); }

inline float dist(const vec3f &a, const vec3f &b)
{ return sqrtf(dist2(a,b)); }

inline bool operator==(const vec3f &a, const vec3f &b)
{ return a.v[0]==b.v[0] && a.v[1]==b.v[1] &&  a.v[2]==b.v[2]; }

inline bool operator!=(const vec3f &a, const vec3f &b)
{ return a.v[0]!=b.v[0] || a.v[1]!=b.v[1] || a.v[2]!=b.v[2]; }

inline vec3f operator-(const vec3f &a)
{ return vec3f(-a.v[0], -a.v[1], -a.v[2]); }

inline vec3f operator+(const vec3f &a, const vec3f &b)
{ return vec3f(a.v[0]+b.v[0], a.v[1]+b.v[1], a.v[2]+b.v[2]); }

inline vec3f operator-(const vec3f &a, const vec3f &b)
{ return vec3f(a.v[0]-b.v[0], a.v[1]-b.v[1], a.v[2]-b.v[2]); }

inline vec3f operator*(const vec3f &a, float scalar)
{ return vec3f(scalar*a.v[0], scalar*a.v[1], scalar*a.v[2]); }


inline vec3f operator*(float scalar, const vec3f &a)
{ return vec3f(scalar*a.v[0], scalar*a.v[1], scalar*a.v[2]); }


inline vec3f operator/(const vec3f &a, float scalar)
{ return vec3f(a.v[0]/scalar, a.v[1]/scalar,  a.v[2]/scalar); }

inline float dot(const vec3f &a, const vec3f &b)
{ return a.v[0]*b.v[0] + a.v[1]*b.v[1] + a.v[2]*b.v[2]; }

inline vec3f cross(const vec3f &a, const vec3f &b)
{ return vec3f(a.v[1]*b.v[2]-b.v[1]*a.v[2], b.v[0]*a.v[2]-a.v[0]*b.v[2], a.v[0]*b.v[1]-b.v[0]*a.v[1]); }

inline void normalize(vec3f &a)
{ a/=mag(a); }

inline vec3f normalized(const vec3f &a)
{ return a/mag(a); }

inline std::ostream &operator<<(std::ostream &out, const vec3f &a)
{ return out<<a.v[0]<<' '<<a.v[1] << ' ' << a.v[2]; }

inline std::istream &operator>>(std::istream &in, vec3f &a)
{ return in>>a.v[0]>>a.v[1]>>a.v[2]; }

}

#endif
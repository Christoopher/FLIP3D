#ifndef _VECTOR3_H_
#define _VECTOR3_H_

#include <iostream>
#include <cmath>
#include "util.h"

struct vec3f
{
	float v[3];
	vec3f(float v_) { v[0] = v[1] = v[2] = 0; }
	vec3f(float x, float y, float z) {v[0] = x; v[1] = y; v[2] = z; }
	vec3f(vec3f & source) { v[0] = source.v[0]; v[1] = source.v[1]; v[2] = source.v[2]; }
	float & operator[](int index) {return v[index]; }
	const vec3f & operator[](int index ) const { return v[index]; }

	vec3f operator+=(const vec3f & rhs)
	{
		v[0] += rhs.v[0];
		v[1] += rhs.v[1];
		v[2] += rhs.v[2];
		return *this;
	}

	vec3f operator-=(const vec3f & rhs)
	{
		v[0] -= rhs.v[0];
		v[1] -= rhs.v[1];
		v[2] -= rhs.v[2];
		return *this;
	}

	vec3f operator*=(const float & scalar)
	{
		v[0] *= scalar;
		v[1] *= scalar;
		v[2] *= scalar;
		return *this;
	}

	vec3f operator/=(const float & scalar)
	{
		v[0] /= scalar;
		v[1] /= scalar;
		v[2] /= scalar;
		return *this;
	}
}; 


inline float mag2(const vec3f &a)
{ return a.v[0]*a.v[0] + a.v[1]*a.v[1] + a.v[2]*a.v[2]; }

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
	



#endif
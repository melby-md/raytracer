#ifndef _VECTOR_HPP
#define _VECTOR_HPP
#include <math.h>

struct Vec3 {
	double x, y, z;
};

static inline Vec3 operator+(Vec3 a, Vec3 b)
{
	return {a.x + b.x, a.y + b.y, a.z + b.z};
}

static inline Vec3 operator-(Vec3 a, Vec3 b)
{
	return {a.x - b.x, a.y - b.y, a.z - b.z};
}

static inline Vec3 operator-(Vec3 v)
{
	return {-v.x, -v.y, -v.z};
}

static inline Vec3 operator*(Vec3 a, Vec3 b)
{
	return {a.x*b.x, a.y*b.y, a.z*b.z};
}

static inline Vec3 operator*(Vec3 v, double s)
{
	return {v.x*s, v.y*s, v.z*s};
}

static inline Vec3 operator*(double s, Vec3 v)
{
	return v*s;
}

static inline Vec3 operator/(Vec3 v, double s)
{
	return {v.x/s, v.y/s, v.z/s};
}

static inline Vec3 operator/(double s, Vec3 v)
{
	return v/s;
}

static inline Vec3 &operator+=(Vec3 &a, Vec3 b)
{
	a = a + b;
	return a;
}

static inline Vec3 &operator-=(Vec3 &a, Vec3 b)
{
	a = a - b;
	return a;
}

static inline Vec3 &operator*=(Vec3 &v, double s)
{
	v = v * s;
	return v;
}

static inline Vec3 &operator*=(Vec3 &a, Vec3 b)
{
	a = a * b;
	return a;
}

static inline Vec3 &operator/=(Vec3 &v, double s)
{
	v = v / s;
	return v;
}

static inline double Dot(Vec3 a, Vec3 b)
{
	return a.x*b.x + a.y*b.y + a.z*b.z;
}

static inline double Length2(Vec3 v)
{
	return Dot(v, v);
}

static inline double Length(Vec3 v)
{
	return sqrt(Length2(v));
}

static inline Vec3 Normalize(Vec3 v)
{
	return v / Length(v);
}

static inline Vec3 Cross(Vec3 a, Vec3 b)
{
	return {
		a.y*b.z - a.z*b.y,
		a.z*b.x - a.x*b.z,
		a.x*b.y - a.y*b.x
	};
}

static inline Vec3 Reflect(Vec3 v, Vec3 n)
{
	return v - 2 * Dot(v, n) * n;
}

static inline Vec3 Refract(Vec3 v, Vec3 n, double eta)
{
	double k = 1 - eta * eta * (1 - Dot(n, v) * Dot(n, v));
	if (k < 0)
		return {};
	else
		return eta * v - (eta * Dot(n, v) + sqrt(k)) * n;
}

#endif

#ifndef _MATHLIB_HPP
#define _MATHLIB_HPP
#include <math.h>
// Why must we suffer Bjarne?

struct Vec2 {
	double x, y;
};

struct Vec3 {
	double x, y, z;
};

static inline Vec3 operator+(Vec3 a, Vec3 b)
{
	return Vec3{a.x + b.x, a.y + b.y, a.z + b.z};
}

static inline Vec3 operator+(Vec3 v, double s)
{
	return v + Vec3{s, s, s};
}

static inline Vec3 operator+(double s, Vec3 v)
{
	return v + s;
}

static inline Vec3 operator-(Vec3 a, Vec3 b)
{
	return Vec3{a.x - b.x, a.y - b.y, a.z - b.z};
}

static inline Vec3 operator-(Vec3 v)
{
	return Vec3{-v.x, -v.y, -v.z};
}

static inline Vec3 operator-(Vec3 v, double s)
{
	return v - Vec3{s, s, s};
}

static inline Vec3 operator-(double s, Vec3 v)
{
	return Vec3{s, s, s} - v;
}

static inline Vec3 operator*(Vec3 a, Vec3 b)
{
	return Vec3{a.x*b.x, a.y*b.y, a.z*b.z};
}

static inline Vec3 operator*(Vec3 v, double s)
{
	return Vec3{v.x*s, v.y*s, v.z*s};
}

static inline Vec3 operator*(double s, Vec3 v)
{
	return v*s;
}

static inline Vec3 operator/(Vec3 a, Vec3 b)
{
	return Vec3{a.x/b.x, a.y/b.y, a.z/b.z};
}

static inline Vec3 operator/(Vec3 v, double s)
{
	return Vec3{v.x/s, v.y/s, v.z/s};
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

static inline Vec3 &operator+=(Vec3 &v, double s)
{
	v = v + s;
	return v;
}

static inline Vec3 &operator-=(Vec3 &a, Vec3 b)
{
	a = a - b;
	return a;
}

static inline Vec3 &operator-=(Vec3 &v, double s)
{
	v = v - s;
	return v;
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

static inline Vec3 &operator/=(Vec3 &a, Vec3 b)
{
	a = a / b;
	return a;
}

static inline Vec3 Lerp(Vec3 v0, Vec3 v1, double t)
{
	return (1 - t) * v0 + t * v1;
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
	return Vec3{
		a.y*b.z - a.z*b.y,
		a.z*b.x - a.x*b.z,
		a.x*b.y - a.y*b.x
	};
}

static inline Vec3 Pow(Vec3 a, Vec3 b)
{
	return Vec3{pow(a.x, b.x), pow(a.y, b.y), pow(a.z, b.z)};
}

static inline Vec3 Exp(Vec3 v)
{
	return Vec3{exp(v.x), exp(v.y), exp(v.z)};
}

static inline Vec3 Reflect(Vec3 v, Vec3 n)
{
	return v - 2 * Dot(v, n) * n;
}

// 3x3 column major
struct Mat3 {
	Vec3 i, j, k;
};

static inline Vec3 operator*(Mat3 m, Vec3 v)
{
	return v.x * m.i + v.y * m.j + v.z * m.k;
}

static inline Mat3 OrthoNormalBasis(Vec3 z_axis)
{
	Vec3 a = (fabs(z_axis.x) > .9) ? Vec3{0, 1, 0} : Vec3{1, 0, 0};
	Vec3 y_axis = Normalize(Cross(z_axis, a));
	Vec3 x_axis = Cross(z_axis, y_axis);

	return Mat3{x_axis, y_axis, z_axis};
}

static inline Mat3 Transpose(Mat3 m)
{
	return Mat3{
		{m.i.x, m.j.x, m.k.x},
		{m.i.y, m.j.y, m.k.y},
		{m.i.z, m.j.z, m.k.z}
	};
}

#endif

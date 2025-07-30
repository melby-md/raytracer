// Why must we suffer Bjarne?
#include <math.h>

constexpr float PI = 3.14159265f;
constexpr float INV_PI = 0.31830988f;

struct Vec2 {
	float x, y;
};

struct Vec3 {
	union {float x; float r;};
	union {float y; float g;};
	union {float z; float b;};
};

static inline float Clamp(float x, float mn, float mx)
{
	return fminf(mx, fmaxf(mn, x));
}

static inline Vec3 Clamp(Vec3 x, float mn, float mx)
{
	return {
		Clamp(x.x, mn, mx),
		Clamp(x.y, mn, mx),
		Clamp(x.z, mn, mx)
	};
}

static inline Vec3 V3(float s)
{
	return {s, s, s};
}

static inline Vec3 operator+(Vec3 a, Vec3 b)
{
	return Vec3{a.x + b.x, a.y + b.y, a.z + b.z};
}

static inline Vec3 operator+(Vec3 v, float s)
{
	return v + Vec3{s, s, s};
}

static inline Vec3 operator+(float s, Vec3 v)
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

static inline Vec3 operator-(Vec3 v, float s)
{
	return v - Vec3{s, s, s};
}

static inline Vec3 operator-(float s, Vec3 v)
{
	return Vec3{s, s, s} - v;
}

static inline Vec3 operator*(Vec3 a, Vec3 b)
{
	return Vec3{a.x*b.x, a.y*b.y, a.z*b.z};
}

static inline Vec3 operator*(Vec3 v, float s)
{
	return Vec3{v.x*s, v.y*s, v.z*s};
}

static inline Vec3 operator*(float s, Vec3 v)
{
	return v*s;
}

static inline Vec3 operator/(Vec3 a, Vec3 b)
{
	return Vec3{a.x/b.x, a.y/b.y, a.z/b.z};
}

static inline Vec3 operator/(Vec3 v, float s)
{
	return Vec3{v.x/s, v.y/s, v.z/s};
}

static inline Vec3 operator/(float s, Vec3 v)
{
	return v/s;
}

static inline Vec3 &operator+=(Vec3 &a, Vec3 b)
{
	a = a + b;
	return a;
}

static inline Vec3 &operator+=(Vec3 &v, float s)
{
	v = v + s;
	return v;
}

static inline Vec3 &operator-=(Vec3 &a, Vec3 b)
{
	a = a - b;
	return a;
}

static inline Vec3 &operator-=(Vec3 &v, float s)
{
	v = v - s;
	return v;
}

static inline Vec3 &operator*=(Vec3 &v, float s)
{
	v = v * s;
	return v;
}

static inline Vec3 &operator*=(Vec3 &a, Vec3 b)
{
	a = a * b;
	return a;
}

static inline Vec3 &operator/=(Vec3 &v, float s)
{
	v = v / s;
	return v;
}

static inline Vec3 &operator/=(Vec3 &a, Vec3 b)
{
	a = a / b;
	return a;
}

static inline Vec3 Lerp(Vec3 v0, Vec3 v1, float t)
{
	return v0 + (v1 - v0) * t;
}

static inline float Dot(Vec3 a, Vec3 b)
{
	return a.x*b.x + a.y*b.y + a.z*b.z;
}

static inline float Length2(Vec3 v)
{
	return Dot(v, v);
}

static inline float Length(Vec3 v)
{
	return sqrtf(Length2(v));
}

static inline Vec3 Normalize(Vec3 v)
{
	float length = Length(v);
	Assert(length > 0);
	return v / length;
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
	return Vec3{powf(a.x, b.x), powf(a.y, b.y), powf(a.z, b.z)};
}

static inline Vec3 Exp(Vec3 v)
{
	return Vec3{expf(v.x), expf(v.y), expf(v.z)};
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
	Vec3 a = (fabsf(z_axis.x) > .9f) ? Vec3{0, 1, 0} : Vec3{1, 0, 0};
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

#define CountOf(x) (sizeof(x)/sizeof((x)[0]))

#define Panic(...) \
	do { \
		Error(__VA_ARGS__); \
		exit(1); \
	} while (0)

#define xstr(x) str(x)
#define str(x) #x

#define _log(level, ...) Log(level ":" __FILE__ ":" xstr(__LINE__) ": " __VA_ARGS__)
#define Error(...) _log("ERROR", __VA_ARGS__)
#define Trace(...) _log("TRACE", __VA_ARGS__)
#define Log(...) fprintf(stderr, __VA_ARGS__)

#ifdef OPT
#  define Assert(c)
#elif __GNUC__
#  define Assert(c) if (!(c)) __builtin_trap()
#elif _MSC_VER
#  define Assert(c) if (!(c)) __debugbreak()
#else
#  undef NDEBUG
#  include <assert.h>
#  define Assert(c) assert(c)
#endif

typedef uintptr_t uptr;
typedef ptrdiff_t isize;
typedef size_t    usize;

typedef uint64_t  u64;
typedef uint32_t  u32;
typedef uint16_t  u16;
typedef uint8_t   u8;

typedef int64_t   i64;
typedef int32_t   i32;
typedef int16_t   i16;
typedef int8_t    i8;

typedef unsigned char byte;

struct String {
	char *data;
	isize length;
};

#define S(s) (String{(char *)(s), sizeof(s)-1})

struct Vec2 {
	float x, y;
};

struct Vec3 {
	union {float x, r;};
	union {float y, g;};
	union {float z, b;};

	float operator[](isize pos)
	{
		Assert(pos >= 0 && pos < 3);
		return ((float *)this)[pos];
	}
};

// 3x3 column major
struct Mat3 {
	Vec3 i, j, k;
};

constexpr float PI = 3.14159265f;

static inline float Min(float a, float b) {
    return a < b ? a : b;
}

static inline float Max(float a, float b) {
    return a > b ? a : b;
}

static inline Vec3 Min(Vec3 a, Vec3 b)
{
	return {
		Min(a.x, b.x),
		Min(a.y, b.y),
		Min(a.z, b.z)
	};
}

static inline Vec3 Max(Vec3 a, Vec3 b)
{
	return {
		Max(a.x, b.x),
		Max(a.y, b.y),
		Max(a.z, b.z)
	};
}

static inline float Clamp(float x, float mn, float mx)
{
	return Min(mx, Max(mn, x));
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
	return {a.x + b.x, a.y + b.y, a.z + b.z};
}

static inline Vec3 operator+(Vec3 v, float s)
{
	return {v.x + s, v.y + s, v.z + s};
}

static inline Vec3 operator+(float s, Vec3 v)
{
	return v + s;
}

static inline Vec3 operator-(Vec3 a, Vec3 b)
{
	return {a.x - b.x, a.y - b.y, a.z - b.z};
}

static inline Vec3 operator-(Vec3 v)
{
	return {-v.x, -v.y, -v.z};
}

static inline Vec3 operator-(Vec3 v, float s)
{
	return {v.x - s, v.y - s, v.z - s};
}

static inline Vec3 operator-(float s, Vec3 v)
{
	return {s - v.x, s - v.y, s - v.z};
}

static inline Vec3 operator*(Vec3 a, Vec3 b)
{
	return {a.x*b.x, a.y*b.y, a.z*b.z};
}

static inline Vec3 operator*(Vec3 v, float s)
{
	return {v.x*s, v.y*s, v.z*s};
}

static inline Vec3 operator*(float s, Vec3 v)
{
	return v*s;
}

static inline Vec3 operator/(Vec3 a, Vec3 b)
{
	return {a.x/b.x, a.y/b.y, a.z/b.z};
}

static inline Vec3 operator/(Vec3 v, float s)
{
	return {v.x/s, v.y/s, v.z/s};
}

static inline Vec3 operator/(float s, Vec3 v)
{
	return {s/v.x, s/v.y, s/v.z};
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

static inline Vec3 operator*(Mat3 m, Vec3 v)
{
	return v.x * m.i + v.y * m.j + v.z * m.k;
}

static inline Mat3 OrthoNormalBasis(Vec3 z_axis)
{
	Vec3 a = (fabsf(z_axis.x) > .9f) ? Vec3{0, 1, 0} : Vec3{1, 0, 0};
	Vec3 y_axis = Normalize(Cross(z_axis, a));
	Vec3 x_axis = Cross(z_axis, y_axis);

	return {x_axis, y_axis, z_axis};
}

static inline Mat3 Transpose(Mat3 m)
{
	return {
		m.i.x, m.j.x, m.k.x,
		m.i.y, m.j.y, m.k.y,
		m.i.z, m.j.z, m.k.z
	};
}

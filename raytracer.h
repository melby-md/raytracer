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

#ifdef RELEASE
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

constexpr float PI = 3.14159265f;
constexpr float INV_PI = 0.31830988f;

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
};

// 3x3 column major
struct Mat3 {
	Vec3 i, j, k;
};

enum ObjectType {
	OBJ_TRIANGLE,
	OBJ_SPHERE
};

struct Light {
	Vec3 color;
	i32 object_idx;
};

struct Triangle {
	Vec3 v0, v1, v2;
	Vec3 n0, n1, n2;
};

struct Sphere {
	Vec3 center;
	float radius;
};

struct Object {
	ObjectType type;
	i32 material_idx;
	i32 light_idx;
	union {
		Triangle triangle;
		Sphere sphere;
	};
};

struct Material {
	Vec3 color;
	float roughness;
	float ior;
	float metallic;
};

struct Sample {
	Vec3 bsdf;
	float pdf;
	Vec3 l;
};

struct Ray {
	Vec3 origin, direction;
};

struct Hit {
	Vec3 point;
	Vec3 normal;
	i32 obj_idx;
};

constexpr i32 MAX_OBJECTS = 64;
constexpr i32 MAX_MATERIALS = 64;
constexpr i32 MAX_LIGHTS = 64;

struct Scene {
	Vec3 camera;
	Vec3 up;
	Vec3 look_at;
	float fov;
	float defocus_angle;
	float exposure;
	i16 width, height;
	i16 samples;

	int object_count;
	Object objects[MAX_OBJECTS];

	int light_count;
	Light lights[MAX_LIGHTS];

	int material_count;
	Material materials[MAX_MATERIALS];

	Vec3 sun_color;
};

// bsdf.cpp
Vec3   BSDF(Vec3, Vec3, Material *);
float  BSDFPDF(Vec3, Vec3, Material *);
Sample SampleBSDF(Vec3, Material *, u32 *);

// main.cpp
float Rand(u32 *);
u32   Rand(u32 *, u32, u32);

// parser.cpp
void LoadScene(Scene *, const char *);

// Math
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
		m.i.x, m.j.x, m.k.x,
		m.i.y, m.j.y, m.k.y,
		m.i.z, m.j.z, m.k.z
	};
}

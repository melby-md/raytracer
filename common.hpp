// Odds and ends...
#include <stddef.h> // ptrdiff_t
#include <stdint.h> // *int*_t
#include <stdio.h>

#define CountOf(x) (sizeof(x)/sizeof((x)[0]))

#define Panic(...) \
	do { \
		Error(__VA_ARGS__); \
		exit(1); \
	} while (0)

#define xstr(x) str(x)
#define str(x) #x

#define _log(level, ...) fprintf(stderr, level ":" __FILE__ ":" xstr(__LINE__) ": " __VA_ARGS__)
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

u32 _rand(u32 *rng_state)
{
	// Xorshift32
	u32 x = *rng_state;
	x ^= x << 13;
	x ^= x >> 17;
	x ^= x << 5;
	*rng_state = x;

	return x;
}

float Rand(u32 *rng_state)
{
	return (float)(_rand(rng_state) >> 8) / (1<<24);
}

u32 Rand(u32 *rng_state, u32 min, u32 max)
{
	return _rand(rng_state) % (max - min + 1) + min;
}

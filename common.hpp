// Odds and ends...
#ifndef _COMMON_HPP
#define _COMMON_HPP
#include <stddef.h> // ptrdiff_t
#include <stdint.h> // *int*_t
#include <stdio.h>

// macros

#define sizeof(x) ((iz)sizeof(x))
#define Countof(x) (sizeof(x)/sizeof((x)[0]))

#define Max(a, b) ((a) > (b) ? (a) : (b))
#define Min(a, b) ((a) < (b) ? (a) : (b))

#define Clamp(x, mn, mx) Min(Max(x, mn), mx)

#define Panic(msg) \
	do { \
		Error("%s", msg); \
		exit(1); \
	} while (0)

#define xstr(x) str(x)
#define str(x) #x

#define _log(level, ...) fprintf(stderr, level ":" __FILE__ ":" xstr(__LINE__) ": " __VA_ARGS__)
#define Error(...) _log("ERROR", __VA_ARGS__)
#define Debug(...) _log("DEBUG", __VA_ARGS__)
#define Log(...) fprintf(stderr, __VA_ARGS__)

#if __GNUC__
#  define Break() __builtin_unreachable()
#elif _MSC_VER
#  define Break() __debugbreak()
#else
#  define Break() (*(volatile int *)0 = 0)
#endif

#define Assert(c) do if (!(c)) Break(); while(0)

// types

typedef uintptr_t uptr;
typedef ptrdiff_t iz;
typedef size_t    uz;

typedef uint64_t  u64;
typedef uint32_t  u32;
typedef uint16_t  u16;
typedef uint8_t   u8;

typedef int64_t   i64;
typedef int32_t   i32;
typedef int16_t   i16;
typedef int8_t    i8;

typedef unsigned char byte;

double Rand(u32 *rng_state)
{
	// Xorshift32
	u32 x = *rng_state;
	x ^= x << 13;
	x ^= x >> 17;
	x ^= x << 5;
	*rng_state = x;

	return (double)x / UINT32_MAX;
}

#endif

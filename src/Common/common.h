#ifndef COMMON_H
#define COMMON_H

#include <stdint.h>

typedef uint8_t  u8;
typedef uint16_t u16;
typedef uint32_t u32;
typedef uint64_t u64;

typedef int8_t  s8;
typedef int16_t s16;
typedef int32_t s32;
typedef int64_t s64;

typedef unsigned char byte;

typedef int bool;

#define true 1
#define false 0

#define ArrayCount(Array) (sizeof(Array) / sizeof(Array[0]))

#define Assert(cond) do { if(!__builtin_expect((cond), 1)) { __builtin_trap(); } } while(0)

#endif // COMMON_H

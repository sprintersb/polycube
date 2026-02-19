#ifndef ENDIAN_H
#define ENDIAN_H

#if !defined ENDIAN_AGNOSTIC
#if __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__
#define ENDIAN_LITTLE
#elif __BYTE_ORDER__ == __ORDER_BIG_ENDIAN__
#define ENDIAN_BIG
#else
#define ENDIAN_AGNOSTIC
#endif // __BYTE_ORDER__
#endif // ! ENDIAN_AGNOSTIC

#endif // ENDIAN_H

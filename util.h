// -*- c++ -*-
#ifndef UTIL_H
#define UTIL_H

#include <cstdint>

namespace gjl
{
    template<typename T>
    constexpr inline T factorial (int n)
    {
        T fac = n > 0 ? n : 1;
        while (--n > 1)
            fac *= n;
        return fac;
    }

    constexpr inline int hyperoctahedral_order (int dim)
    {
        return factorial<int> (dim) * (1 << dim);
    }

    template <typename T>
    constexpr inline T bitswap (T t, int a, int b)
    {
        const T maska = (T) 1 << a;
        const T maskb = (T) 1 << b;
        const T mask = maska ^ maskb;
        const T tm = t & mask;
        return tm == 0 || tm == mask ? t : t ^ mask;
    }

    // Only used with non-zero arguments.
    template <typename T>
    constexpr inline int count_trailing_zeros (T t);
    template<> constexpr inline int count_trailing_zeros (short t) { return __builtin_ctz (t); }
    template<> constexpr inline int count_trailing_zeros (unsigned short t) { return __builtin_ctz (t); }
    template<> constexpr inline int count_trailing_zeros (int t) { return __builtin_ctz (t); }
    template<> constexpr inline int count_trailing_zeros (unsigned t) { return __builtin_ctz (t); }
    template<> constexpr inline int count_trailing_zeros (long t) { return __builtin_ctzl (t); }
    template<> constexpr inline int count_trailing_zeros (unsigned long t) { return __builtin_ctzl (t); }
    template<> constexpr inline int count_trailing_zeros (long long t) { return __builtin_ctzll (t); }
    template<> constexpr inline int count_trailing_zeros (unsigned long long t) { return __builtin_ctzll (t); }
#if __SIZEOF_LONG_LONG__ != 32 && defined DIM && DIM > 8
    template<> constexpr inline int count_trailing_zeros (unsigned __int128 t)
    {
        const uint64_t hi = (uint64_t) (t >> 64);
        const uint64_t lo = (uint64_t) t;
        return lo != 0
            ? count_trailing_zeros (lo)
            : 64 + count_trailing_zeros (hi);
    }
#endif

// Only used with non-zero arguments.
    template <typename T>
    constexpr inline int count_leading_zeros (T t);
    template<> constexpr inline int count_leading_zeros (short t) { return __builtin_clz (t); }
    template<> constexpr inline int count_leading_zeros (unsigned short t) { return __builtin_clz (t); }
    template<> constexpr inline int count_leading_zeros (int t) { return __builtin_clz (t); }
    template<> constexpr inline int count_leading_zeros (unsigned t) { return __builtin_clz (t); }
    template<> constexpr inline int count_leading_zeros (long t) { return __builtin_clzl (t); }
    template<> constexpr inline int count_leading_zeros (unsigned long t) { return __builtin_clzl (t); }
    template<> constexpr inline int count_leading_zeros (long long t) { return __builtin_clzll (t); }
    template<> constexpr inline int count_leading_zeros (unsigned long long t) { return __builtin_clzll (t); }
#if __SIZEOF_LONG_LONG__ != 32 && defined DIM && DIM > 8
    template<> constexpr inline int count_leading_zeros (unsigned __int128 t)
    {
        const uint64_t hi = (uint64_t) (t >> 64);
        const uint64_t lo = (uint64_t) t;
        return hi != 0
            ? count_leading_zeros (hi)
            : 64 + count_leading_zeros (lo);
    }
#endif

    template <typename T>
    constexpr inline int popcount (T t);
    template<> constexpr inline int popcount (int t) { return __builtin_popcount (t); }
    template<> constexpr inline int popcount (unsigned t) { return __builtin_popcount (t); }
    template<> constexpr inline int popcount (long t) { return __builtin_popcountl (t); }
    template<> constexpr inline int popcount (unsigned long t) { return __builtin_popcountl (t); }
    template<> constexpr inline int popcount (long long t) { return __builtin_popcountll (t); }
    template<> constexpr inline int popcount (unsigned long long t) { return __builtin_popcountll (t); }

    // When t0 has a popcount of 2, then return the LSB and MSB bit positions.
    template <typename T>
    constexpr inline bool popcount2_bits (T t0, int &lsb, int &msb)
    {
        const T t1 = t0 & (t0 - 1);
        const T t2 = t1 & (t1 - 1);
        if (! t1 || t2)
            return false;
        lsb = gjl::count_trailing_zeros (t0);
        msb = gjl::count_trailing_zeros (t1);
        return true;
    }
}; // gjl

#endif // UTIL_H


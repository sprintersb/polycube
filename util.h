// -*- c++ -*-
#ifndef UTIL_H
#define UTIL_H

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

template <typename T>
constexpr inline int count_trailing_zeros (T t);

template<> constexpr inline int count_trailing_zeros (int t) { return __builtin_ctz (t); }
template<> constexpr inline int count_trailing_zeros (unsigned t) { return __builtin_ctz (t); }
template<> constexpr inline int count_trailing_zeros (long t) { return __builtin_ctzl (t); }
template<> constexpr inline int count_trailing_zeros (unsigned long t) { return __builtin_ctzl (t); }
template<> constexpr inline int count_trailing_zeros (long long t) { return __builtin_ctzll (t); }
template<> constexpr inline int count_trailing_zeros (unsigned long long t) { return __builtin_ctzll (t); }

#endif // UTIL_H


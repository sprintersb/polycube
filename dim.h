// -*- c++ -*-
#ifndef DIM_H
#define DIM_H

#include <algorithm> // std::min
#include <iostream>
#include <array>
// C
#include <cstdint>
#include <cassert>
// Own
#include "hash.h"

// Make CUBES_ARRAY / DimArray the default for performance.
#if defined CUBES_ARRAY && defined CUBES_VECT
#error pick one of -DCUBES_ARRAY or -DCUBES_VECT
#elif !defined CUBES_ARRAY && !defined CUBES_VECT
#define CUBES_ARRAY
#elif defined CUBES_VECT
#undef CELLS
#endif

#if defined CUBES_ARRAY && !defined CELLS
#error CUBES_ARRAY requires to define -DCELLS=<value>
#endif

struct DimIterator;

struct Dim
{
    using value_t = int8_t;
#if DIM == 2
    using vector_t = value_t __attribute__((vector_size(2)));
    using int_t = uint16_t;
#elif DIM > 2 && DIM <= 4
    using vector_t = value_t __attribute__((vector_size(4)));
    using int_t = uint32_t;
#elif DIM > 4 && DIM <= 8
    using vector_t = value_t __attribute__((vector_size(8)));
    using int_t = uint64_t;
#elif DIM > 8 && DIM <= 16
    using vector_t = value_t __attribute__((vector_size(16)));
    using int_t = unsigned __int128;
#else
#error DIM=?
#endif

    vector_t v = (vector_t) (int_t) 0;
    static inline constexpr vector_t all0 = (vector_t) (int_t) 0;
    Dim () {}
    Dim (std::initializer_list<value_t> &&il)
    {
        assert (il.size () == 0 || il.size () == DIM);
        v = Dim::all0;
        int j = 0;
        for (auto i : il)
            v[j++] = i;
    }
    Dim (const vector_t &v) : v(v) {}

    static Dim all (int w)
    {
        Dim d{};
        for (int i = 0; i < d.size (); ++i)
            d.v[i] = w;
        return d;
    }
    int_t ival () const
    {
        return (int_t) v;
    }
    static inline int size ()
    {
        return DIM;
    }
    void set (int i, int val)
    {
#if defined CUBES_ARRAY
        // v[i] = Gives warning with CUBES_ARRAY.
        vector_t w(v);
        w[i] = (Dim::value_t) val;
        v = w;
#else
        v[i] = (Dim::value_t) val;
#endif
    }
    int operator [] (int i) const
    {
        return v[i];
    }
    // Shift-invariant comparison, so (int_t) <=> (int_t) won't work.
    // Benefit is that Cubes.cells don't change their order when shifted.
    int cmp (Dim d) const
    {
        for (int i = 0; i < size (); ++i)
            if (v[i] != d.v[i])
                return v[i] - d.v[i];
        return 0;
    }
    bool operator == (Dim d) const { return (int_t) v == (int_t) d.v; }
    bool operator != (Dim d) const { return (int_t) v != (int_t) d.v; }
    bool operator <= (Dim d) const { return cmp (d) <= 0; }
    bool operator >= (Dim d) const { return cmp (d) >= 0; }
    bool operator <  (Dim d) const { return cmp (d) <  0; }
    bool operator >  (Dim d) const { return cmp (d) >  0; }
    void operator += (Dim d) { v += d.v; }
    void operator -= (Dim d) { v -= d.v; }
    void operator *= (int i) { v *= (value_t) i; }
    Dim operator + (Dim d) const { return Dim (v + d.v); }
    Dim operator - (Dim d) const { return Dim (v - d.v); }
    Dim operator * (int i) const { return Dim (v * (value_t) i); }
    int operator % (Dim d) const { return v[0] * d[1] - v[1] * d[0]; }
    int operator * (Dim d) const
    {
        int s = 0;
        for (int i = 0; i < size (); ++i)
            s += v[i] * d[i];
        return s;
    }
    Dim rot (int i /* Left in units of 90 deg */) const
    {
        Dim d (*this);
        i = (4 + (i % 4)) % 4;
        while (i-- > 0)
            d = Dim { (value_t) -d.v[1], d.v[0] };
        return d;
    }
    void min (Dim d)
    {
        for (int i = 0; i < size (); ++i)
            v[i] = std::min (v[i], d.v[i]);
    }
    void max (Dim d)
    {
        for (int i = 0; i < size (); ++i)
            v[i] = std::max (v[i], d.v[i]);
    }
    Hasher::type hash () const
    {
        return Hasher() (ival ());
    }
    struct Hash
    {
        Hasher::type operator () (Dim d) const
        {
            return d.hash ();
        }
    };
    // Manhattan
    int abs () const
    {
        int man = std::abs (v[0]);
        for (int i = 1; i < size (); ++i)
            man += std::abs (v[i]);
        return man;
    }
    // Manhattan distance
    int dist (Dim r) const
    {
        return (*this - r).abs ();
    }
    static Dim rand (int lo, int hi)
    {
        const int mod = hi - lo + 1;
        Dim d;
        for (int i = 0; i < size (); ++i)
            d.set (i, lo + std::rand () % mod);
        return d;
    }
    // Cube diagonal starting at vertex id, pointing to the opposite vertex.
    static Dim diag (int id)
    {
        Dim d;
        for (int i = 0; i < size (); ++i)
            d.set (i, (id & (1 << i)) ? -1 : 1);
        return d;
    }
    // The distance to the line p0 + <v>, multiplied by DIM.
    int dist (Dim p0, Dim v) const
    {
        const Dim d = *this - p0;
        const Dim r = d * size() - v * (d * v);
        return r.abs ();
    }
    static Dim rand (Dim range)
    {
        Dim d;
        for (int i = 0; i < size (); ++i)
            d.set (i, std::rand () % range[i]);
        return d;
    }

    DimIterator begin () const;
    DimIterator end () const;
    static const Dim Min;
    static const Dim Max;
    friend std::ostream& operator << (std::ostream&, Dim);
};

inline const Dim Dim::Min = Dim::all (INT8_MIN);
inline const Dim Dim::Max = Dim::all (INT8_MAX);

struct Box
{
    Dim lo, hi;
    bool contains (Dim d)
    {
        for (int i = 0; i < d.size (); ++i)
            if (d.v[i] < lo.v[i] || d.v[i] > hi.v[i])
                return false;
        return true;
    }
    Dim vertex (int id) const
    {
        Dim d;
        for (int i = 0; i < lo.size (); ++i)
            d.set (i, (id & (1 << i)) ? hi[i] : lo[i]);
        return d;
    }
};


class DimIterator
{
    using Corona0 = std::array<Dim, 2 * DIM>;
    static inline Corona0 get_corona0 ()
    {
        Corona0 c;
        for (int i = 0; i < 2 * DIM; ++i)
            c[i].v[i / 2] = i % 2 == 0 ? 1 : -1;
        return c;
    }
    static inline const Corona0 corona0 = DimIterator::get_corona0 ();
    int pos;

    DimIterator (int pos) : pos(pos) {}
    friend Dim;
public:
    bool operator != (DimIterator di) const
    {
        return pos != di.pos;
    }
    void operator ++ ()
    {
        ++ pos;
    }
    Dim operator * () const
    {
        return DimIterator::corona0[pos];
    }
};

inline DimIterator Dim::begin () const { return DimIterator (0); }
inline DimIterator Dim::end ()   const { return DimIterator (2 * DIM); }

inline std::ostream& operator << (std::ostream &ost, Dim d)
{
    char comma = '<';
    for (int i = 0; i < d.size (); ++i)
    {
        ost << comma << (int) d[i];
        comma = ',';
    }
    return ost << '>';
}

#ifdef CELLS
class DimArray
{
    std::array<Dim, 1 + CELLS> a_;
public:
    Dim operator [] (int i) const
    {
        return a_[i];
    }
    int size () const
    {
        int sz = (int) (Dim::int_t) a_[CELLS].v;
        assert (sz >= 0 && sz <= CELLS);
        return sz;
    }
    void clear ()
    {
        set_size (0);
    }
    struct Iterator
    {
        friend DimArray;
        void operator ++ () { ++ptr; };
        bool operator == (const Iterator &i) const { return ptr == i.ptr; }
        bool operator != (const Iterator &i) const { return ptr != i.ptr; }
        Dim& operator * () const { return *ptr; }
    private:
        Dim *ptr;
        Iterator () = delete;
        Iterator (Dim *ptr) : ptr(ptr) {}
    };
    struct CIterator
    {
        friend DimArray;
        void operator ++ () { ++ptr; };
        bool operator == (const CIterator &i) const { return ptr == i.ptr; }
        bool operator != (const CIterator &i) const { return ptr != i.ptr; }
        Dim operator * () const { return *ptr; }
    private:
        const Dim *ptr;
        CIterator () = delete;
        CIterator (const Dim *ptr) : ptr(ptr) {}
    };
    Iterator begin () { return Iterator (&a_[0]); }
    Iterator end ()   { return Iterator (&a_[size ()]); }
    CIterator begin () const { return CIterator (&a_[0]); }
    CIterator end ()   const { return CIterator (&a_[size ()]); }
    CIterator cbegin () const { return CIterator (&a_[0]); }
    CIterator cend ()   const { return CIterator (&a_[size ()]); }
    void insert (const Iterator &it, Dim d)
    {
        const int pos = (int) (& (*it) - & a_[0]);
        assert (pos >= 0 && pos <= size ());
        for (int i = size (); i > pos; --i)
            a_[i] = a_[i - 1];
        a_[pos] = d;
        set_size (1 + size ());
    }
private:
    void set_size (int sz)
    {
        assert (sz >= 0 && sz <= CELLS);
        a_[CELLS].v = (Dim::vector_t) (Dim::int_t) sz;
    }
};
#endif // CELLS

#endif // DIM_H

// -*- c++ -*-
#ifndef DIM_H
#define DIM_H

#include "config-dim.h" // auto-generated
#ifndef DIM
#error DIM is not defined
#endif

#include <algorithm> // std::min
#include <iostream>
#include <array>
#include <iterator>
// C
#include <cstdint>
#include <cstring>
// Own
#include "hash.h"
#include "iterator-wrap.h"
#include "util.h"
#include "endian.h"
#include "debug.h"

// Make CUBES_ARRAY / DimArray the default for performance.
#if defined CUBES_ARRAY && defined CUBES_VECT
#error pick one of -DCUBES_ARRAY or -DCUBES_VECT
#elif !defined CUBES_ARRAY && !defined CUBES_VECT
#define CUBES_ARRAY
#endif

#if defined CUBES_ARRAY && !defined CELLS
#error CUBES_ARRAY requires to define -DCELLS=<value>
#endif

template<typename T>
using VertexTraits = std::array<T, 1 << DIM>;
struct DimIterator;

struct Dim
{
    using value_t = int8_t;
#if DIM == 2
    using vector_t = [[gnu::vector_size (2)]] value_t;
    using int_t = uint16_t;
#elif DIM > 2 && DIM <= 4
    using vector_t = [[gnu::vector_size (4)]] value_t;
    using int_t = uint32_t;
#elif DIM > 4 && DIM <= 8
    using vector_t = [[gnu::vector_size (8)]] value_t;
    using int_t = uint64_t;
#elif DIM > 8 && DIM <= 16
    using vector_t = [[gnu::vector_size (16)]] value_t;
    using int_t = unsigned __int128;
#else
#error DIM=?
#endif

    vector_t v = all0;
    static inline constexpr vector_t all0 = (vector_t) (int_t) 0;
    static const vector_t mask;
    static const VertexTraits<Dim> diagonals;
    static const VertexTraits<Dim> vertex_masks;
    Dim () {}
    Dim (const std::initializer_list<value_t> &il)
    {
        Assert (il.size () == 0 || il.size () == DIM,
                "bad initializer length %zd", il.size ());
        v = Dim::all0;
        int j = 0;
        for (auto i : il)
            v[j++] = i;
    }
    Dim (const vector_t &v) : v(v) {}

    static inline Dim all (int w)
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
        v[i] = (Dim::value_t) val;
    }
    int operator [] (int i) const
    {
        return v[i];
    }
    static inline vector_t combine (vector_t cond, vector_t v1, vector_t v0)
    {
        return (v1 & cond) | (v0 & ~cond);
    }
    static inline Dim combine (Dim cond, Dim v1, Dim v0)
    {
        return Dim (Dim::combine (cond.v, v1.v, v0.v));
    }
    // Shift-invariant comparison, so (int_t) <=> (int_t) won't work.
    // Benefit is that Cubes.cells don't change their order when shifted.
    int cmp (Dim d) const
    {
        d -= *this;
#if DIM == 2
        return d[0] ? d[0] : d[1];
#elif defined ENDIAN_AGNOSTIC
        for (int i = 0; i < size () - 1; ++i)
            if (d[i])
                return d[i];
        return d[size () - 1];
#else // Known endianess and DIM >= 3.
        int_t i = d.ival ();
        if (i == 0)
            return 0;
        // Avoid a loop under the assumption that CTZ / CLZ is fast.
#ifdef ENDIAN_LITTLE
        const int lsb = gjl::count_trailing_zeros (i);
        i >>= lsb & ~7;
#elif defined ENDIAN_BIG
        const int msb = gjl::count_leading_zeros (i);
        i <<= msb & ~7;
#else
#error unknown endianess
#endif // Little or Big.
        return ((vector_t) i)[0];
#endif // Known endianess.
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
#if DIM == 2
    void min (Dim d)
    {
        v[0] = std::min (v[0], d.v[0]);
        v[1] = std::min (v[1], d.v[1]);
    }
    void max (Dim d)
    {
        v[0] = std::max (v[0], d.v[0]);
        v[1] = std::max (v[1], d.v[1]);
    }
#else
    void min (Dim d)
    {
        v = Dim::combine (v < d.v, v, d.v);
    }
    void max (Dim d)
    {
        v = Dim::combine (v > d.v, v, d.v);
    }
#endif // DIM
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
    int trace () const
    {
        int tr = v[0];
        for (int i = 1; i < size (); ++i)
            tr += v[i];
        return tr;
    }
    // Manhattan
    int abs () const
    {
#if DIM == 2
        return std::abs (v[0]) + std::abs (v[1]);
#elif DIM == 3
        // Whether this is profitable for 3D depends on the machine.
        return std::abs (v[0]) + std::abs (v[1]) + std::abs (v[2]);
#else
        return Dim (Dim::combine (v < Dim::all0, -v, v)).trace ();
#endif
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
    static inline Dim diag (int id)
    {
        return diagonals[id];
    }
    static inline Dim vertex_mask (int id)
    {
        return vertex_masks[id];
    }
private:
    // Cube diagonal starting at vertex id, pointing to the opposite vertex.
    static Dim make_diag (int id)
    {
        Dim d;
        for (int i = 0; i < size (); ++i)
            d.set (i, (id & (1 << i)) ? -1 : 1);
        return d;
    }
    static Dim make_mask (int id)
    {
        Dim d;
        for (int i = 0; i < size (); ++i)
            d.set (i, (id & (1 << i)) ? 0 : 0xff);
        return d;
    }
public:
    // The distance to the line p0 + <v>, multiplied by DIM = |v|^2.
    int dist (Dim p0, Dim v) const
    {
        const Dim d = *this - p0;
        const Dim r = d * size() - v * (d * v);
        return r.abs ();
    }
    int pseudo_dist (Dim p0, Dim) const
    {
        const Dim d = *this - p0;
        return d * d;
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
inline const Dim::vector_t Dim::mask = Dim::make_mask ((1 << DIM) - 1).v;
inline const VertexTraits<Dim> Dim::diagonals =
    [] ()
    {
        VertexTraits<Dim> vt;
        for (size_t id = 0; id < vt.size (); ++id)
            vt[id] = Dim::make_diag (id);
        return vt;
    } ();
inline const VertexTraits<Dim> Dim::vertex_masks =
    [] ()
    {
        VertexTraits<Dim> vt;
        for (size_t id = 0; id < vt.size (); ++id)
            vt[id] = Dim::make_mask (id);
        return vt;
    } ();


struct Box
{
    Dim lo, hi;
    bool operator == (const Box &r) const
    {
        return lo == r.lo && hi == r.hi;
    }
    bool contains (Dim d) const
    {
        for (int i = 0; i < d.size (); ++i)
            if (d.v[i] < lo.v[i] || d.v[i] > hi.v[i])
                return false;
        return true;
    }
    Dim vertex (int id) const
    {
        return Dim::combine (Dim::vertex_mask (id), lo, hi);
    }
};


// Iterate over the neighbors of 0.  This means that iteration only depends
// on DIM and not on the actual location of the producer.
class DimIterator
{
    using Corona0 = std::array<Dim, 2 * DIM>;
    static inline const Corona0 corona0 = []()
    {
        Corona0 c;
        for (int i = 0; i < 2 * DIM; ++i)
            c[i].v[i / 2] = i % 2 == 0 ? 1 : -1;
        return c;
    } ();
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

#ifdef CUBES_ARRAY
// We will be at the brink of memory exhaustion, so we need to save
// each possible byte on storage and therefore use a managed std::array.
// When space and time don't matter, you can still use a std::vector by
// re-compiling with -DCUBES_VECT.
class DimArray
{
    std::array<Dim, 1 + CELLS> a_;
public:
    using value_type = Dim;
    Dim* data ()
    {
        return &a_[0];
    }
    const Dim* data () const
    {
        return &a_[0];
    }
    Dim operator [] (int i) const
    {
        return a_[i];
    }
    int size () const
    {
        const int sz = (int) a_[CELLS].ival ();
        Assert (sz >= 0 && sz <= CELLS, "get bad size %d", sz);
        return sz;
    }
    void clear ()
    {
        set_size (0);
    }

    using iterator       = IteratorWrap<DimArray, Dim*, Dim&>;
    using const_iterator = IteratorWrap<DimArray, const Dim*, Dim>;

    auto begin () { return iterator (&a_[0]); }
    auto end ()   { return iterator (&a_[size ()]); }
    auto cbegin () const { return const_iterator (&a_[0]); }
    auto cend ()   const { return const_iterator (&a_[size ()]); }
    auto begin () const { return cbegin (); }
    auto end ()   const { return cend (); }
    void insert (const iterator &it, Dim d)
    {
        const int pos = (int) (& (*it) - & a_[0]);
        Assert (pos >= 0 && pos <= size (), "bad insert position %d", pos);
        const auto n_bytes = (const uint8_t*) &*end() - (const uint8_t*) &*it;
        if (n_bytes < 16)
            for (int i = size (); i > pos; --i)
                a_[i] = a_[i - 1];
        else
            std::memmove (&*(it + 1), &*it, n_bytes);
        *it = d;
        set_size (1 + size ());
    }
private:
    void set_size (int sz)
    {
        Assert (sz >= 0 && sz <= CELLS, "set bad size %d", sz);
        a_[CELLS].v = (Dim::vector_t) (Dim::int_t) sz;
    }
};

template<>
struct std::iterator_traits<DimArray::iterator>
{
    using iterator_category = forward_iterator_tag;
    using value_type = Dim; // DimArray::iterator::value_type;
    using difference_type = ptrdiff_t;
};

#endif // CUBES_ARRAY

#endif // DIM_H

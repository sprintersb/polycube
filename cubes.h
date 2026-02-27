// -*- c++ -*-
#ifndef CUBES_H
#define CUBES_H

#include <string>
#include <iostream>
#include <utility> // std::swap
#include <algorithm> // std::sort
#include <atomic>
// Containers
#include <optional>
#include <array>
#include <vector>
#include <set>
// Own
#include "hash.h"
#include "dim.h"
#include "iterator-wrap.h"
#include "corona.h"
#include "debug.h"

// The maximal possible length (in Manhattan metric) of the diagonal
// of a Cubes' bounding box.  This is used during canonicalization.
#if !defined MAX_DIAGONAL_LENGTH
#ifdef CELLS
#define MAX_DIAGONAL_LENGTH CELLS
#else
#warning using MAX_DIAGONAL_LENGTH = 80
#define MAX_DIAGONAL_LENGTH 80
#endif
#endif // MAX_DIAGONAL_LENGTH

// For sizes from SORT_THRESHOLD on, garbled cells are sorted with std::sort.
// For sizes below, sorting is achieved by re-building the cells.
#ifndef SORT_THRESHOLD
#define SORT_THRESHOLD 11
#endif

// For sizes from BINARY_ADD_THRESHOLD on the insert position is determined
// by a binary search.  For sizes below a linear search is used.
#ifndef BINARY_ADD_THRESHOLD
#define BINARY_ADD_THRESHOLD 10
#endif

// The relative cost of a vertex canonicalization compared to one
// step of traversing the hyperoctahedral group.  For stats only.
#define VERTEX_CANONICALIZATION_COST 4

enum
{
    STAT_FAIL, STAT_SUCC, STAT_COST_FAIL, STAT_COST_SUCC,
    STAT_sentinel
};

inline std::array<std::atomic<int64_t>, STAT_sentinel> stat;

struct Dist;
using VertexValues = std::array<int, 1 << DIM>;
using DistPointers = std::array<Dist*, 1 << DIM>;

// Cubes is a vector / array since that is most memory friendly.
// Normalized such that bounding box lower is all 0's.
struct Cubes
{
private:
#ifdef CUBES_ARRAY
    using container_type = DimArray;
#else
    using container_type = std::vector<Dim>;
#endif
    container_type cells;

public:
    static constexpr bool can_canonicalize = 2 <= DIM && DIM <= 5;
    static inline bool take_stat;
    // Mirror symmetry along the specified dimension.
    using Symmetry = std::optional<int>;

    // The id of a vertex in [0, 2^DIM).
    using Vertex = std::optional<int>;

    struct Context;
    Cubes () {}
    Cubes (const std::initializer_list<Dim> &il)
    {
        for (Dim d : il)
            add (d);
    }
    Cubes (const Cubes &dad, Dim d)
    {
        cells = dad.cells;
        add (d);
    }

    int size () const
    {
        return cells.size ();
    }
    int cmp (const Cubes &r) const
    {
        if (size () != r.size ())
            return size () - r.size ();
        auto p2 = r.cells.begin ();
        for (Dim d : cells)
            if (const int i = d.cmp (*p2); i != 0)
                return i;
            else
                ++p2;
        return 0;
    }
    bool operator == (const Cubes &r) const
    {
        if (size () != r.size ())
            return false;
        auto p2 = r.cells.begin ();
        for (Dim d : cells)
            if (d != *p2)
                return false;
            else
                ++p2;
        return true;
    }
    bool operator < (const Cubes &r) const
    {
        return this->cmp (r) < 0;
    }
    Hasher::type hash () const
    {
        Hasher::type h = 0;
        for (Dim d : cells)
            h = Hasher::add (h, d.ival());
        return h;
    }
    bool contains (Dim d) const
    {
        for (Dim c : cells)
            if (int i = c.cmp (d); i == 0)
                return true;
            else if (i > 0)
                break;
        return false;
    }
    void add_search_linear (Dim d)
    {
        for (auto it = cells.begin (); ; ++it)
            if (int i; it == cells.end () || (i = d.cmp (*it)) < 0)
            {
                cells.insert (it, d);
                break;
            }
            else
                Assert (i != 0, "assume we always increase cells");
    }
    void add_search_binary (Dim d)
    {
        int l = 0;
        for (int r = size (); r > l; )
        {
            const int m = (l + r) >> 1;
            if (const int i = d.cmp (cells[m]); i > 0)
                l = m + 1;
            else if (i < 0)
                r = m;
            else
                Assert (i != 0, "assume we always increase cells");
        }
        cells.insert (cells.begin() + l, d);
    }
    void add (Dim d)
    {
        if (size () >= BINARY_ADD_THRESHOLD)
            add_search_binary (d);
        else
            add_search_linear (d);
        // Normalize again.  Max one component of d is negative.
        for (int j = 0; j < d.size (); ++j)
            if (d.v[j] < 0)
            {
                shift (j, -d.v[j]);
                break;
            }
    }
    // Notice that shift() is compatible with Dim.cmp(), hence
    // order will be retained.
    void shift (int i, int off)
    {
        for (Dim &d : cells)
            d.set (i, d[i] + off);
    }
    int max_coord (int i) const
    {
        int8_t m = INT8_MIN;
        for (Dim d : cells)
            m = std::max<int> (m, d[i]);
        return m;
    }
public:
    Cubes rotated (int i, int j) const
    {
        Cubes c (*this);
        return c.rotate (i, j);
    }
    Cubes& rotate (int i, int j)
    {
        const int xm = max_coord (i);
        if (size () < SORT_THRESHOLD)
        {
            Cubes c;
            for (Dim d : cells)
            {
                const int y = xm - d[i];
                d.set (i, d[j]);
                d.set (j, y);
                c.add (d);
            }
            cells = std::move (c.cells);
        }
        else
        {
            for (Dim &d : cells)
            {
                const int y = xm - d[i];
                d.set (i, d[j]);
                d.set (j, y);
            }
            sort ();
        }
        return *this;
    }
    Cubes mirrored (int i) const
    {
        const int m = max_coord (i);
        Cubes c;
        for (Dim d : cells)
        {
            d.set (i, m - d[i]);
            c.add (d);
        }
        return c;
    }
    int squeeze (Box &bbox)
    {
        int l = 0, r = DIM - 1;
        for (bool swapped = false; ; swapped = true)
        {
            while (r > 0 && bbox.hi[r] == 0)  --r;
            while (l < r && bbox.hi[l] != 0)  ++l;
            if (l >= r)
            {
                if (swapped)
                    sort ();
                return 1 + r;
            }
            swap (l, r);
            std::swap (bbox.hi.v[l], bbox.hi.v[r]);
        }
    }
    // This is usually not needed since add() adds Dim's in order.
    // It must be run after bulk changes like in rotate() or squeeze().
    void sort ()
    {
        std::sort (cells.begin (), cells.end ());
    }

    using iterator =
        IteratorWrap<container_type, container_type::iterator, Dim&>;
    using const_iterator =
        IteratorWrap<container_type, container_type::const_iterator, Dim>;
    auto begin () { return iterator (cells.begin ()); }
    auto end ()   { return iterator (cells.end ()); }
    auto cbegin () const { return const_iterator (cells.cbegin ()); }
    auto cend () const   { return const_iterator (cells.cend ()); }
    auto begin () const { return cbegin (); }
    auto end () const   { return cend (); }
private:
    // Flip all the dimensions as specified in mask.
    // !!! This will trash the ordering !!!
    void flip (int mask, const Box &bbox)
    {
        Assert (bbox.lo == Dim::all0, "non-aligned Cubes");
        for (Dim &d : *this)
            for (int b = 0; b < DIM; ++b)
                if (mask & (1 << b))
                    d.set (b, bbox.hi[b] - d[b]);
    }
    // !!! This will trash the ordering !!!
    void flip (int mask)
    {
        flip (mask, bounding_box ());
    }
    // Swap the specified dimensions.  !!! This will trash the ordering !!!
    void swap (int a, int b)
    {
        for (Dim &d : *this)
        {
            const int ia = d[a];
            const int ib = d[b];
            d.set (a, ib);
            d.set (b, ia);
        }
    }
    bool matches_flipped (int dim) const
    {
        return matches_flipped (dim, max_coord (dim));
    }
    bool matches_flipped (int dim, const Box &bbox) const
    {
        return matches_flipped (dim, bbox.hi[dim]);
    }
    bool matches_flipped (int dim, int hi_dim) const
    {
        for (Dim d : *this)
        {
            d.set (dim, hi_dim - d[dim]);
            if (! contains (d))
                return false;
        }
        return true;
    }
public:
    void canonicalize ();
    void canonicalize (bool &symmetric);
    int multiplicity () const;
    Cubes canonical () const
    {
        Cubes c (*this);
        c.canonicalize ();
        return c;
    }
    bool is_canonical () const
    {
        return *this == canonical ();
    }
private:
    template<typename T>
    struct CongruentsAspect {};
    template<typename T>
    T congruents (int dim, bool same_parity) const;
    bool canonical_vertex (Context&) const;
    bool maybe_canonicalize_vertices (Context&, bool same_parity = false);
    void canonicalize_vertices (Context&, bool);
    Symmetry find_symmetry (const DistPointers&, const Box&, int) const;
public:
    Box bounding_box () const
    {
        Box box { size() ? cells[0] : Dim::Max, size() ? cells[0] : Dim::Min };
        for (int i = 1; i < size (); ++i)
        {
            box.lo.min (cells[i]);
            box.hi.max (cells[i]);
        }
        return box;
    }
    Corona corona () const
    {
        Corona cora;
        for (Dim d : *this)
            for (Dim delta : d)
                if (! Cubes::contains (d + delta))
                    cora.add (d + delta);
        return cora;
    }
    std::string ascii (char c = '*') const;
    friend std::ostream& operator << (std::ostream&, const Cubes&);
};

struct Cubes::Context
{
    int dim;
    Box bbox;
    Vertex vertex;
    VertexValues vvs;
    Symmetry symmetry;
};

#endif // CUBES_H

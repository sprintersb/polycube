// -*- c++ -*-
#ifndef CUBES_H
#define CUBES_H

#if defined CUBES_ARRAY && !defined CELLS
#error CUBES_ARRAY CELLS=?
#endif

#include <string>
#include <iostream>
#include <vector>
#include <set>
#include <utility> // std::swap
// C
#include <cstring>
#include <cstdlib>
#include <omp.h>
// Own
#include "hash.h"
#include "dim.h"
#include "util.h"

struct CubesIterator;
struct CubesMuIterator;
struct Pad;
struct Dist;
using VertexValues = std::array<int, 1 << DIM>;

// Cubes is a vector / array since that is most memory friendly.
struct Cubes
{
private:
    // Normalized such that bounding box lower is all 0's.
#ifdef CUBES_ARRAY
    using cell_iterator = DimArray::CIterator;
    using cell_muiterator = DimArray::Iterator;
    DimArray cells;
#else
    using cell_iterator = std::vector<Dim>::const_iterator;
    using cell_muiterator = std::vector<Dim>::iterator;
    std::vector<Dim> cells;
#endif

public:
    Cubes () {}
    Cubes (const Cubes *dad, Dim d)
    {
        if (dad)
            cells = dad->cells;
        add (d);
    }
    int size () const
    {
        return cells.size ();
    }
    int cmp (const Cubes &c) const
    {
        auto &&p2 = c.cells.begin ();
        auto &&e2 = c.cells.end ();
        for (Dim d : cells)
        {
            if (p2 == e2)
                return 1;
            const int i = d.cmp (*p2);
            if (i)
                return i;
            ++p2;
        }
        return p2 == e2 ? 0 : -1;
    }
    bool operator == (const Cubes &r) const
    {
        return this->cmp (r) == 0;
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
    void add (Dim d)
    {
        for (auto it = cells.begin (); ; ++it)
            if (int i; it == cells.end () || (i = d.cmp (*it)) < 0)
            {
                cells.insert (it, d);
                break;
            }
            else if (i == 0)
                assert (i != 0 && "Assume we always increase cells");
        // Normalize again.  Max one component of d is negative.
        for (int j = 0; j < d.size (); ++j)
            if (d.v[j] < 0)
            {
                shift (j, -d.v[j]);
                break;
            }
    }
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
        const int xm = max_coord (i);
        Cubes c;
        for (Dim d : cells)
        {
            int y = xm - d[i];
            d.set (i, d[j]);
            d.set (j, y);
            c.add (d);
        }
        return c;
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
private:
    // Flip all the dimensions as specified in mask.
    // !!! This will trash the ordering !!!
    void flip (int mask, const Box&);
    // Swap the specified dimensions.  !!! This will trash the ordering !!!
    void swap (int, int);
    bool matches_flipped (int dim, const Box&) const;
public:
    Pad congruents () const;
    void canonicalize ();
    int multiplicity () const;
    Cubes canonical () const
    {
        Cubes c (*this);
        c.canonicalize ();
        return c;
    }
    int maybe_canonical_vertex (VertexValues&, const Box&, int&) const;
    bool canonicalizable_vertices (int &symmetry) const;
    bool maybe_canonicalize_vertices ();
    int maybe_symmetry (const std::array<Dist*, 1 << DIM>&, const Box&) const;
    CubesMuIterator begin ();
    CubesMuIterator end ();
    CubesIterator begin () const;
    CubesIterator end () const;
    CubesIterator cbegin () const;
    CubesIterator cend () const;
    friend CubesIterator;
    friend CubesMuIterator;

    // Common tail.
public:
    Box bounding_box () const;
    std::string ascii (char c = '*') const;
    friend std::ostream& operator << (std::ostream&, const Cubes&);
};

struct CubesIterator
{
    Cubes::cell_iterator it;
    void operator ++ () { ++it; }
    bool operator == (const CubesIterator &r) const { return it == r.it; }
    bool operator != (const CubesIterator &r) const { return it != r.it; }
    Dim operator * () const { return *it; }
};

struct CubesMuIterator
{
    Cubes::cell_muiterator it;
    void operator ++ () { ++it; }
    bool operator == (const CubesMuIterator &r) const { return it == r.it; }
    bool operator != (const CubesMuIterator &r) const { return it != r.it; }
    Dim& operator * () const { return *it; }
};

inline CubesIterator Cubes::cbegin () const
{
    return CubesIterator { cells.cbegin () };
}

inline CubesIterator Cubes::cend () const
{
    return CubesIterator { cells.cend () };
}

inline CubesIterator Cubes::begin () const { return cbegin (); }
inline CubesIterator Cubes::end () const   { return cend (); }

inline CubesMuIterator Cubes::begin ()
{
    return CubesMuIterator { cells.begin () };
}

inline CubesMuIterator Cubes::end ()
{
    return CubesMuIterator { cells.end () };
}

inline void Cubes::flip (int mask, const Box &bbox)
{
    assert (bbox.lo == Dim::all0);
    for (Dim &d : *this)
        for (int b = 0; b < DIM; ++b)
            if (mask & (1 << b))
                d.set (b, bbox.hi[b] - d[b]);
}

inline bool Cubes::matches_flipped (int dim, const Box &bbox) const
{
    for (Dim d : *this)
    {
        d.set (dim, bbox.hi[dim] - d[dim]);
        if (! contains (d))
            return false;
    }
    return true;
}

inline void Cubes::swap (int a, int b)
{
    for (Dim &d : *this)
    {
        const int ia = d[a];
        const int ib = d[b];
        d.set (a, ib);
        d.set (b, ia);
    }
}

struct Pad : std::set<Cubes> {};

Pad Cubes::congruents () const
{
    Cubes c (*this);
    Pad pad;
    const char *S = "";
#if DIM == 2
    S = "111" "0" "111";
#elif DIM == 3
    S = ("12121 3 12121 3 12121 1 12121 0"
         "12121 3 12121 3 12121 1 12121");
#elif DIM == 4
#define M "11121113111311141112111"
    S = (M "5" M "5" M "6" M "5" M "6" M "5" M "5" M "0"
         M "5" M "5" M "6" M "5" M "6" M "5" M "5" M);
#undef M
#else
    assert (0 && "todo: canonicalize in DIM");
#endif
    pad.insert (c);
    for (auto s = S; *s; ++s)
        switch (*s - '0')
        {
#define PADD(FF) pad.insert (c = c.FF); break
            default: continue;
            case 0: PADD (mirrored (0));
            case 1: PADD (rotated (0, 1));
#if DIM == 3
            case 2: PADD (rotated (1, 2));
            case 3: PADD (rotated (2, 0));
#elif DIM == 4
            case 2: PADD (rotated (0, 2));
            case 3: PADD (rotated (1, 2));
            case 4: PADD (rotated (2, 1));
            case 5: PADD (rotated (3, 0));
            case 6: PADD (rotated (3, 1));
#endif
#undef PADD
        } // switch
    return pad;
}

std::array<std::atomic<int64_t>, 3> stat;

void Cubes::canonicalize ()
{
    bool success = maybe_canonicalize_vertices ();
    stat[success] += 1;
    if (! success)
        cells = std::move (congruents().begin()->cells);
}

int Cubes::multiplicity () const
{
    int symmetry;
    return canonicalizable_vertices (symmetry)
        ? hyperoctahedral_order (DIM) >> (symmetry >= 0)
        : (int) congruents().size ();
}


// Keeping track of how many cublis have a specific distance to
// a vertex of the bounding box.  The lengths of the diagonals
// of our bounding boxes are all below 80 (in Manhattan metric).
using DistBase = std::array<int16_t, 80>;
struct Dist : DistBase
{
    int id = -1;
    int value = -1;
    Dist ()
    {
        std::memset (DistBase::data(), 0,
                     DistBase::size() * sizeof (DistBase::value_type));
    }
    int size () const
    {
        return size_;
    }
    value_type& at (int i)
    {
        if (i >= size_)
            resize (1 + i);
        return DistBase::operator [] (i);
    }
    value_type at (int i) const
    {
        return i >= size_ ? 0 : DistBase::operator [] (i);
    }
    value_type& operator [] (int i) { return at (i); }
    value_type  operator [] (int i) const { return at (i); }
    int cmp (const Dist &r) const
    {
        if (size () != r.size ())
            return size () - r.size ();
        for (int i = 0; i < size (); ++i)
            if (at (i) != r.at (i))
                return at (i) - r.at (i);
        return 0;
    }
    bool neighbors_uniquely_sortable (const VertexValues &vvs) const
    {
        unsigned bits = 0;
        int n_bits = 0;
        for (int j = 0; j < DIM; ++j)
        {
            const unsigned mask = 1u << vvs[id ^ (1 << j)];
            n_bits += ! (bits & mask);
            bits |= mask;
        }
        return n_bits == DIM;
    }
private:
    int size_ = 0;
    void resize (int i)
    {
        assert (i < (int) DistBase::size ());
        size_ = i;
    }
    iterator begin () = delete;
    iterator end   () = delete;
    const_iterator begin () const = delete;
    const_iterator end   () const = delete;
    const_iterator cbegin () const = delete;
    const_iterator cend   () const = delete;
};

std::ostream& operator << (std::ostream &ost, const Dist &c)
{
    ost << "Dist." << c.id << ":  \t";
    for (int i = 0; i < c.size (); ++i)
        ost << " " << + c.at (i);
    return ost << "\n";
}

inline bool Cubes::canonicalizable_vertices (int &symmetry) const
{
    VertexValues vvs;
    return maybe_canonical_vertex (vvs, bounding_box (), symmetry) >= 0;
}

// Return dim when there is respective simple mirror symmetry,
// otherwise return -1.
int Cubes::maybe_symmetry (const std::array<Dist*, 1 << DIM> &pc,
                           const Box &bbox) const
{
    const int max_value = pc[(1 << DIM) - 1]->value;
    if (max_value != (1 << (DIM - 1)) - 1)
        return -1;
    // # valuations are exacly half of HOh's symmetries, which is a
    // hallmark for mirror symmetry.
    const int symmetry = pc[0]->id ^ pc[1]->id;
    // For now, only consider simple reflections, i.e. symmetry has popcount 1.
    if (symmetry & (symmetry - 1))
        return -1;
    for (int i = 2; i < 1 << DIM; i += 2)
        if ((pc[i]->id ^ pc[i + 1]->id) != symmetry)
            return -1;
    // All looks good until now, but we need a final proof.
    const int dim = count_trailing_zeros (symmetry);
    if (! matches_flipped (dim, bbox))
        return -1;
    return dim;
}

// Return the id of a canonical vertex in [0, 2^DIM), or -1 if not found.
// Set symmetry to a dimension of mirror symmetry, or -1 if none such.
int Cubes::maybe_canonical_vertex (VertexValues &vvs, const Box &bbox,
                                   int &symmetry) const
{
    symmetry = -1;
    // For each cubli, record its man distances to any of the bounding vertices.
    std::array<Dist, 1 << DIM> vdist;
    for (int id = 0; id < 1 << DIM; ++id)
    {
        // The id of a vertex is the canonical integer in [0, 2^DIM).
        Dist &vd = vdist[id];
        vd.id = id;
        const Dim diag = Dim::diag (id);
        const Dim ecke = bbox.vertex (id);
        for (Dim d : *this)
            // Get the distance to bounding vertex id's cube diagonal.
            // This is better than just adding 1 for a cubli of its distance.
            vd[d.dist (ecke)] += 1 + 13 * d.dist (ecke, diag);
    }

    // Sort vdist[].
    std::array<Dist*, 1 << DIM> pc;
    for (int id = 0; id < 1 << DIM; ++id)
        pc[id] = & vdist[id];
    // Prefer qsort over std::sort since qsort can use spaceship (cmp).
    std::qsort (&pc, pc.size (), sizeof (Dist*),
                [](const void *va, const void *vb) -> int
                {
                    const Dist *pa = *(const Dist *const *) va;
                    const Dist *pb = *(const Dist *const *) vb;
                    return pa->cmp (*pb);
                });
    // Add valuations in increasing order, assigning same valuation
    // if it matches the predecessor's.
    pc[0]->value = 0;
    for (size_t b = 1; b < 1 << DIM; ++b)
        pc[b]->value = pc[b - 1]->value + (pc[b]->cmp (*pc[b - 1]) > 0);
    const int max_value = pc[(1 << DIM) - 1]->value;
    // We are after the vertex -> valuation.  Transfer them from auto to output.
    for (size_t id = 0; id < 1 << DIM; ++id)
        vvs[id] = vdist[id].value;
    // Shortcut: When all vertex dists are different, we always succeed.
    if (max_value == (1 << DIM) - 1)
        // Canonical: Use the vertex with the smallest valuation.
        return pc[0]->id;

    if ((symmetry = maybe_symmetry (pc, bbox)) >= 0)
        return pc[0]->id;

    // Now most likely we see a Cubes with some symmetry, but in up
    // to 10% of cases we succeed by looking a bit more closely.
    for (int b = 0; b < 1 << DIM; ++b)
    {
        bool good = true;
        good &= b == 0              || pc[b]->value > pc[b - 1]->value;
        good &= b == (1 << DIM) - 1 || pc[b]->value < pc[b + 1]->value;
        if (good && pc[b]->neighbors_uniquely_sortable (vvs))
            // Canonical: The vertex with the smallest unique valuation.
            return pc[b]->id;
    }
    return -1;
}

bool Cubes::maybe_canonicalize_vertices ()
{
    const Box bbox = bounding_box ();
    assert (bbox.lo == Dim::all0);
    VertexValues vvs;
    int symmetry;
    const int id = maybe_canonical_vertex (vvs, bbox, symmetry);
    if (id < 0)
        return false;

    assert (vvs.size () == 1 << DIM);
    Cubes p2 (*this);
    // Reflect such that vertex id becomes vertex 0.
    if (id != 0)
    {
        p2.flip (id, bbox);
        // Adjust vvs accordingly.
        for (int b = 0; b < 1 << DIM; ++b)
            if (b > (b ^ id))
                std::swap (vvs[b], vvs[b ^ id]);
    }
    // Now sort the DIM vertices adjacent to id.
    uint64_t todo = 0;
    for (int b = 0; b < DIM; ++b)
        todo |= (uint64_t) 1 << vvs[1 << b];
    // Vertex 0 has DIM neighbors, canonicalize them such that higher
    // dimension neighbors get higher vertex values.
    for (int dim = 0; dim < DIM - 1; ++dim, todo &= todo - 1)
    {
        assert (todo);
        const int lsb = count_trailing_zeros (todo);
        // Search the vertex with the value as specified by todo in
        // the vertex -> value list. We only have to look for vertex
        // indices that are a power of 2 (and hence are a neighbors of 0).
        for (int b = 0; b < DIM; ++b)
            if (vvs[1 << b] == lsb
                && b != dim)
            {
                // Make vertex b the next neighbor of vertex 0.
                p2.swap (dim, b);
                // Adjust vvs accordingly.
                for (int i = 0; i < DIM; ++i)
                {
                    const int u = 1 << i;
                    const int v = bitswap (u, dim, b);
                    if (u > v)
                        std::swap (vvs[u], vvs[v]);
                }
            }
    }
    // Now p2 has canonical cublis, but swapping and flipping
    // clobbered cells' order.  Re-construct a proper one.
    cells.clear ();
    for (Dim d : p2)
        add (d);

    return true;
}

inline Box Cubes::bounding_box () const
{
    Box box { Dim::Max, Dim::Min };
    for (Dim d : *this)
    {
        box.lo.min (d);
        box.hi.max (d);
    }
    return box;
}

inline std::string Cubes::ascii (char c) const
{
#if DIM == 2
    auto bbox = bounding_box ();
    std::string str;
    for (Dim::value_t y = bbox.hi[1]; y >= bbox.lo[1]; --y)
    {
        for (Dim::value_t x = bbox.lo[0]; x <= bbox.hi[0]; ++x)
            str += contains (Dim { x, y }) ? c : ' ';
        str += "\n";
    }
    return str;
#else
    (void) c;
    return "Cubes.ascii(DIM != 2)";
#endif
}

inline std::ostream& operator << (std::ostream &ost, const Cubes &c)
{
    ost << "{#" << c.size ();
    for (Dim d : c)
        ost << " " << d;
    return ost << " }";
}


#include "cubes-border.h"

#endif // CUBES_H

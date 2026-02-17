// -*- c++ -*-
#ifndef CUBES_H
#define CUBES_H

#include <string>
#include <iostream>
#include <array>
#include <vector>
#include <set>
#include <utility> // std::swap
#include <algorithm> // std::sort
// C
#include <cstring>
#include <cstdlib>
#include <omp.h>
// Own
#include "hash.h"
#include "dim.h"
#include "util.h"
#include "iterator-wrap.h"

// The maximal possible length (in Manhattan metric) of the diagonal
// of a Cubes' bounding box.  This is used during canonicalization.
#if !defined MAX_DIAGONAL_LENGTH
#ifdef CELLS
#define MAX_DIAGONAL_LENGTH CELLS
#else
#define MAX_DIAGONAL_LENGTH 80
#endif
#endif // MAX_DIAGONAL_LENGTH

// For sizes up to SORT_THRESHOLD we sort by re-building he cells.
// For values above we use std::sort on garbled cells.
#ifndef SORT_THRESHOLD
#define SORT_THRESHOLD 14
#endif

// The relative cost of a vertex canonicalization compared to one
// step of traversing the hyperoctahedral group.  For stats only.
#define VERTEX_CANONICALIZATION_COST 4

enum
{
    STAT_FAIL, STAT_SUCC, STAT_COST_FAIL, STAT_COST_SUCC,
    STAT_sentinel
};
std::array<std::atomic<int64_t>, STAT_sentinel> stat;

struct Pad;
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
        auto p2 = c.cells.begin ();
        auto e2 = c.cells.end ();
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
        Cubes c (*this);
        return c.rotate (i, j);
    }
    Cubes& rotate (int i, int j)
    {
        const int xm = max_coord (i);
        if (size () <= SORT_THRESHOLD)
        {
            Cubes c;
            for (Dim d : cells)
            {
                int y = xm - d[i];
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
                int y = xm - d[i];
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
private:
    // Flip all the dimensions as specified in mask.
    // !!! This will trash the ordering !!!
    void flip (int mask, const Box&);
    // Swap the specified dimensions.  !!! This will trash the ordering !!!
    void swap (int, int);
    bool matches_flipped (int dim, const Box&) const;
public:
    Pad congruents (int dim) const;
    void canonicalize ();
    int multiplicity () const;
    Cubes canonical () const
    {
        Cubes c (*this);
        c.canonicalize ();
        return c;
    }
    int maybe_canonical_vertex (VertexValues&, const Box&, int dim, int&) const;
    bool canonicalizable_vertices (int dim, int &symmetry) const;
    bool maybe_canonicalize_vertices (const Box&, int dim);
    int maybe_symmetry (const DistPointers&, const Box&, int) const;

    using iterator =
        IteratorWrap<container_type, container_type::iterator, Dim&>;
    using const_iterator =
        IteratorWrap<container_type, container_type::const_iterator, Dim>;
    iterator begin () { return iterator (cells.begin ()); }
    iterator end ()   { return iterator (cells.end ()); }
    const_iterator cbegin () const { return const_iterator (cells.cbegin ()); }
    const_iterator cend () const   { return const_iterator (cells.cend ()); }
    const_iterator begin () const { return cbegin (); }
    const_iterator end () const   { return cend (); }
public:
    Box bounding_box () const;
    std::string ascii (char c = '*') const;
    friend std::ostream& operator << (std::ostream&, const Cubes&);
};

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

// Mn is used to encode a traversal of the hyperoctahedral group Oh of
// dimension n.  Each letter stands for an elemment of Oh that yields an
// element of Oh not yet encountered.  There is no known formula for such a
// traversal that works in any dimension, hence use the tight encoding below.
// Plus, this is in the hot path, and we want to avoid dups at any cost.
#define M2  "111"
#define M3  M2 "2" M2 "3" M2 "3" M2 "4" M2 "2" M2
#define M4  M3 "5" M3 "5" M3 "6" M3 "5" M3 "6" M3 "5" M3 "5" M3
#define M5  M4 "8" M4 "9" M4 "8" M4 "9" M4 "8" M4 "8" M4 "7" M4 "7" M4 "7" M4

Pad Cubes::congruents (int dim) const
{
    assert (dim >= 1 && dim <= DIM && "todo: DIM");
    Cubes c (*this);
    Pad pad;
    pad.insert (c);
    const char *S = "?";
    switch (dim)
    {
        case 1: S = "0"; break;
        case 2: S = M2 "0" M2; break;
        case 3: S = M3 "0" M3; break;
        case 4: S = M4 "0" M4; break;
        case 5: S = M5 "0" M5; break;
        default:
            assert (0 && "todo: canonicalize in DIM");
    }
    for (auto s = S; *s; ++s)
        switch (*s - '0')
        {
            default: assert (0 && "bad char");
            case 0: pad.insert (c = c.mirrored (0)); break;
            case 1: pad.insert (c.rotate (0, 1)); break;
#if DIM >= 3
            case 2: pad.insert (c.rotate (0, 2)); break;
            case 3: pad.insert (c.rotate (1, 2)); break;
            case 4: pad.insert (c.rotate (2, 1)); break;
#endif
#if DIM >= 4
            case 5: pad.insert (c.rotate (3, 0)); break;
            case 6: pad.insert (c.rotate (3, 1)); break;
#endif
#if DIM >= 5
            case 7: pad.insert (c.rotate (4, 0)); break;
            case 8: pad.insert (c.rotate (4, 1)); break;
            case 9: pad.insert (c.rotate (3, 4)); break;
#endif
        } // switch
    return pad;
}

// Keeping track of how many cublis have a specific distance to a vertex
// of the bounding box.  The lengths of the diagonals of our bounding boxes
// are all below MAX_DIAGONAL_LENGTH (in Manhattan metric).
// This struct is in the hot path, so we use a managed std::array instead
// of a more dynamic container like std::vector or std::map.
using DistBase = std::array<int16_t, 1 + MAX_DIAGONAL_LENGTH>;
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
    bool neighbors_uniquely_sortable (const VertexValues &vvs, int dim) const
    {
        unsigned bits = 0;
        for (int j = 0; j < dim; ++j)
            if (const unsigned mask = 1u << vvs[id ^ (1 << j)]; bits & mask)
                return false;
            else
                bits |= mask;
        return true;
    }
private:
    int size_ = 0;
    void resize (int i)
    {
        assert (i < (int) DistBase::size () && "see MAX_DIAGONAL_LENGTH");
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

// Return the dimension number when there is respective simple
// mirror symmetry, otherwise return -1.
int Cubes::maybe_symmetry (const DistPointers &pc,
                           const Box &bbox, int dim) const
{
    const int max_value = pc[(1 << dim) - 1]->value;
    if (max_value != (1 << (dim - 1)) - 1)
        return -1;
    // # valuations are exacly half of HOh's symmetries, which is a
    // hallmark for mirror symmetry.
    const int symmetry = pc[0]->id ^ pc[1]->id;
    // For now, only consider simple reflections, i.e. symmetry has popcount 1.
    if (symmetry & (symmetry - 1))
        return -1;
    for (int i = 2; i < 1 << dim; i += 2)
        if ((pc[i]->id ^ pc[i + 1]->id) != symmetry)
            return -1;
    // All looks good until now, but we need a final proof.
    const int d = count_trailing_zeros (symmetry);
    return matches_flipped (d, bbox) ? d : -1;
}

int Cubes::multiplicity () const
{
    int symmetry;
    if (canonicalizable_vertices (DIM, symmetry))
        return hyperoctahedral_order (DIM) >> (symmetry >= 0);
    if (DIM >= 3)
    {
        Cubes c (*this);
        Box bbox = c.bounding_box ();
        const int dim = std::max (1, c.squeeze (bbox));
        if (dim == DIM - 1 && c.canonicalizable_vertices (dim, symmetry))
            return DIM * (hyperoctahedral_order (dim) >> (symmetry >= 0));
    }
    return (int) congruents (DIM).size ();
}

inline bool Cubes::canonicalizable_vertices (int dim, int &symmetry) const
{
    VertexValues vvs;
    return maybe_canonical_vertex (vvs, bounding_box (), dim, symmetry) >= 0;
}

void Cubes::canonicalize ()
{
    Box bbox = bounding_box ();
    assert (bbox.lo == Dim::all0 && "expecting aligned cubes");
    assert (bbox.hi.dist (bbox.lo) < (int) std::tuple_size<DistBase>::value
            && "Cubes diameter exceeds DistBase capacity");
    const int dim = DIM >= 3
        ? std::max (1, squeeze (bbox))
        : DIM;
    const bool success = maybe_canonicalize_vertices (bbox, dim);
    stat[success] += 1;
    stat[2 + success] += success
        ? VERTEX_CANONICALIZATION_COST + (dim != DIM)
        : hyperoctahedral_order (dim) + (dim != DIM);
    if (! success)
        cells = std::move (congruents (dim).begin()->cells);
}

// Return the id of a canonical vertex in [0, 2^DIM), or -1 if not found.
// Set symmetry to a dimension of mirror symmetry, or -1 if none such.
int Cubes::maybe_canonical_vertex (VertexValues &vvs, const Box &bbox,
                                   int dim, int &symmetry) const
{
    symmetry = -1;
    // For each cubli, record its man distances to any of the bounding vertices.
    std::array<Dist, 1 << DIM> vdist;
    for (int id = 0; id < 1 << dim; ++id)
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
    DistPointers pc;
    for (int id = 0; id < 1 << DIM /* sic! to avoid warn */; ++id)
        pc[id] = & vdist[id];
    std::sort (pc.begin (), pc.begin () + (1 << dim),
               [] (const Dist *pa, const Dist *pb) -> bool
               {
                   return pa->cmp (*pb) < 0;
               });
    // Add valuations in increasing order, assigning same valuation
    // if it matches the predecessor's.
    pc[0]->value = 0;
    for (int b = 1; b < 1 << dim; ++b)
        pc[b]->value = pc[b - 1]->value + (pc[b]->cmp (*pc[b - 1]) > 0);
    const int max_value = pc[(1 << dim) - 1]->value;
    // We are after the vertex -> valuation.  Transfer them from auto to output.
    for (int id = 0; id < 1 << dim; ++id)
        vvs[id] = vdist[id].value;
    // Shortcut: When all vertex dists are different, we always succeed.
    if (max_value == (1 << dim) - 1)
        // Canonical: Use the vertex with the smallest valuation.
        return pc[0]->id;

    if ((symmetry = maybe_symmetry (pc, bbox, dim)) >= 0)
        return pc[0]->id;

    // Now most likely we see a Cubes with some symmetry, but in up
    // to 10% of cases we succeed by looking a bit more closely.
    for (int b = 0; b < 1 << dim; ++b)
    {
        bool good = true;
        good &= b == 0              || pc[b]->value > pc[b - 1]->value;
        good &= b == (1 << dim) - 1 || pc[b]->value < pc[b + 1]->value;
        if (good && pc[b]->neighbors_uniquely_sortable (vvs, dim))
            // Canonical: The vertex with the smallest unique valuation.
            return pc[b]->id;
    }
    return -1;
}

bool Cubes::maybe_canonicalize_vertices (const Box &bbox, int dim)
{
    VertexValues vvs;
    int symmetry;
    const int id = maybe_canonical_vertex (vvs, bbox, dim, symmetry);
    if (id < 0)
        return false;

    assert (vvs.size () == 1 << DIM);
    Cubes p2 (*this);
    // Reflect such that vertex id becomes vertex 0.
    if (id != 0)
    {
        p2.flip (id, bbox);
        // Adjust vvs accordingly.
        for (int b = 0; b < 1 << dim; ++b)
            if (b > (b ^ id))
                std::swap (vvs[b], vvs[b ^ id]);
    }
    // Now sort the dim vertices adjacent to id.
    uint64_t todo = 0;
    for (int b = 0; b < dim; ++b)
        todo |= (uint64_t) 1 << vvs[1 << b];
    // Vertex 0 has dim neighbors, canonicalize them such that higher
    // dimension neighbors get higher vertex values.
    for (int d = 0; d < dim - 1; ++d, todo &= todo - 1)
    {
        assert (todo);
        const int lsb = count_trailing_zeros (todo);
        // Search the vertex with the value as specified by todo in
        // the vertex -> value list. We only have to look for vertex
        // indices that are a power of 2 (and hence are a neighbors of 0).
        for (int b = 0; b < dim; ++b)
            if (vvs[1 << b] == lsb
                && b != d)
            {
                // Make vertex b the next neighbor of vertex 0.
                p2.swap (d, b);
                // Adjust vvs accordingly.
                for (int i = 0; i < dim; ++i)
                {
                    const int u = 1 << i;
                    const int v = bitswap (u, d, b);
                    if (u > v)
                        std::swap (vvs[u], vvs[v]);
                }
            }
    }
    // Now p2 has canonical cublis, but swapping and flipping
    // clobbered cells' order.  Re-construct a proper one.
    if (size () <= SORT_THRESHOLD)
    {
        cells.clear ();
        for (Dim d : p2)
            add (d);
    }
    else
    {
        cells = std::move (p2.cells);
        sort ();
    }

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

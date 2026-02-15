// -*- c++ -*-
#ifndef CUBES_VECT_H
#define CUBES_VECT_H

#include <string>
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

#ifndef CUBES_H
#error use #include "cubes.h"
#endif

struct CubesVectIterator;
struct CubesVectMuIterator;
struct CubesRel;
struct Pad;
using VertexValues = std::array<int, 1 << DIM>;

// CubesVect is a vector / array since that is most memory friendly.
struct CubesVect
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
    CubesVect () {}
    CubesVect (const CubesVect *dad, Dim d)
    {
        if (dad)
            cells = dad->cells;
        add (d);
    }
    int size () const
    {
        return cells.size ();
    }
    int cmp (const CubesVect &c) const
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
    bool operator == (const CubesVect &r) const
    {
        return this->cmp (r) == 0;
    }
    bool operator < (const CubesVect &r) const
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
//private:
    friend CubesRel;
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
    CubesVect rotated (int i, int j) const
    {
        const int xm = max_coord (i);
        CubesVect c;
        for (Dim d : cells)
        {
            int y = xm - d[i];
            d.set (i, d[j]);
            d.set (j, y);
            c.add (d);
        }
        return c;
    }
    CubesVect mirrored (int i) const
    {
        const int m = max_coord (i);
        CubesVect c;
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
public:
    Pad congruents () const;
    void canonicalize ();
    int multiplicity () const;
    CubesVect canonical () const
    {
        CubesVect c (*this);
        c.canonicalize ();
        return c;
    }
    int canonical_vertex (VertexValues&, const Box&) const;
    bool canonicalizable_vertices () const;
    bool maybe_canonicalize_vertices ();
    CubesVectMuIterator begin ();
    CubesVectMuIterator end ();
    CubesVectIterator begin () const;
    CubesVectIterator end () const;
    CubesVectIterator cbegin () const;
    CubesVectIterator cend () const;
    friend CubesVectIterator;
    friend CubesVectMuIterator;

    // Common tail.
public:
    Box bounding_box () const;
    std::string ascii (char c = '*') const;
    friend std::ostream& operator << (std::ostream&, const CubesVect&);
};

struct CubesVectIterator
{
    CubesVect::cell_iterator it;
    void operator ++ () { ++it; }
    bool operator == (const CubesVectIterator &r) const { return it == r.it; }
    bool operator != (const CubesVectIterator &r) const { return it != r.it; }
    Dim operator * () const { return *it; }
};

struct CubesVectMuIterator
{
    CubesVect::cell_muiterator it;
    void operator ++ () { ++it; }
    bool operator == (const CubesVectMuIterator &r) const { return it == r.it; }
    bool operator != (const CubesVectMuIterator &r) const { return it != r.it; }
    Dim& operator * () const { return *it; }
};

inline CubesVectIterator CubesVect::cbegin () const
{
    return CubesVectIterator { cells.cbegin () };
}

inline CubesVectIterator CubesVect::cend () const
{
    return CubesVectIterator { cells.cend () };
}

inline CubesVectIterator CubesVect::begin () const { return cbegin (); }
inline CubesVectIterator CubesVect::end () const   { return cend (); }

inline CubesVectMuIterator CubesVect::begin ()
{
    return CubesVectMuIterator { cells.begin () };
}

inline CubesVectMuIterator CubesVect::end ()
{
    return CubesVectMuIterator { cells.end () };
}

inline void CubesVect::flip (int mask, const Box &bb)
{
    assert (bb.lo == Dim::all0);
    for (Dim &d : *this)
        for (int b = 0; b < DIM; ++b)
            if (mask & (1 << b))
                d.set (b, bb.hi[b] - d[b]);
}

inline void CubesVect::swap (int a, int b)
{
    for (Dim &d : *this)
    {
        const int ia = d[a];
        const int ib = d[b];
        d.set (a, ib);
        d.set (b, ia);
    }
}

#if DIM <= 3
struct Pad : std::vector<CubesVect>
{
    bool add (const CubesVect &c)
    {
        for (auto it = begin (); ; ++it)
            if (int i; it == end () || (i = c.cmp (*it)) < 0)
            {
                std::vector<CubesVect>::insert (it, c);
                return true;
            }
            else if (i == 0)
                break;
        return false;
    }
};
#else
struct Pad : std::set<CubesVect>
{
    bool add (const CubesVect &c)
    {
        const auto sz = size ();
        insert (c);
        return sz != size ();
    }
    void reserve (int) const {}
    void resize (int n, const CubesVect &c)
    {
        assert (n == 1);
        clear ();
        add (c);
    }
    const_reference front () const
    {
        return * cbegin ();
    }
};
#endif

Pad CubesVect::congruents () const
{
    CubesVect c (*this);
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
    pad.reserve (hyperoctahedral_order (DIM));
    pad.resize (1, c);
    for (auto s = S; *s; ++s)
        switch (*s - '0')
        {
#define PADD(FF) pad.add (c = c.FF); break
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

void CubesVect::canonicalize ()
{
    bool success = maybe_canonicalize_vertices ();
    stat[success] += 1;
    if (! success)
        cells = std::move (congruents().front().cells);
}

int CubesVect::multiplicity () const
{
    return canonicalizable_vertices ()
        ? hyperoctahedral_order (DIM)
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

inline bool CubesVect::canonicalizable_vertices () const
{
    VertexValues vvs;
    return canonical_vertex (vvs, bounding_box ()) >= 0;
}

// Return the id of a canonica vertex in [0, 2^DIM), or -1 if not found.
// The latter typically occurs when there is some sort of symmetry.
int CubesVect::canonical_vertex (VertexValues &vvs, const Box &bb) const
{
    // For each cubli, record its man distances to any of the bounding vertices.
    std::array<Dist, 1 << DIM> vdist;
    for (int id = 0; id < 1 << DIM; ++id)
    {
        // The id of a vertex is the canonical integer in [0, 2^DIM).
        Dist &vd = vdist[id];
        vd.id = id;
        const Dim diag = Dim::diag (id);
        const Dim ecke = bb.vertex (id);
        for (Dim d : *this)
        {
            // Get the distance to bounding vertex id's cube diagonal.
            // This is better than just adding 1 for a cubli of its distance.
            int dist = 10 + 10 * d.dist (ecke, diag);
            vd[d.dist (ecke)] += 1 + dist * dist;
        }
    }

    // Sorting the VertexDist.
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
    if (1 + max_value == 1 << DIM)
        // Canonical: Use the vertex with the smallest valuation.
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

bool CubesVect::maybe_canonicalize_vertices ()
{
    const Box bb = bounding_box ();
    VertexValues vvs;
    const int id = canonical_vertex (vvs, bb);
    if (id < 0)
        return false;

    assert (vvs.size () == 1 << DIM);
    CubesVect p2 (*this);
    assert (bb.lo == Dim::all0);
    // Reflect such that vertex id becomes vertex 0.
    p2.flip (id, bb);
    // Adjust vvs accordingly.
    for (int b = 0; b < 1 << DIM; ++b)
        if (b > (b ^ id))
            std::swap (vvs[b], vvs[b ^ id]);
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
    // Now p2 is canonical, but swapping and flipping clobbered
    // cells' order.  Re-construct a proper one.
    cells.clear ();
    for (Dim d : p2)
        add (d);

    return true;
}


#undef  CUBES
#define CUBES CubesVect
#include "cubes-common.def"

#endif // CUBES_VECT_H

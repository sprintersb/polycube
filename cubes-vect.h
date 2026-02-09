// -*- c++ -*-
#ifndef CUBES_VECT_H
#define CUBES_VECT_H

// CubesVect is a vect since that is most memory friendly.

#include <string>
#include <vector>
#include <omp.h>
// Own
#include "hash.h"
#include "dim.h"

#ifndef CUBES_H
#error use #include "cubes.h"
#endif

struct CubesVectIterator;
struct CubesRel;

struct CubesVect
{
private:
    // Normalized such that bounding box lower is all 0's.
#ifdef CUBES_ARRAY
    using cell_iterator = DimArray::CIterator;
    DimArray cells;
#else
    using cell_iterator = std::vector<Dim>::const_iterator;
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
    CubesVect rotated () const
    {
        const int xm = max_coord (0);
        CubesVect c;
        for (Dim d : cells)
        {
            int y = xm - d[0];
            d.set (0, d[1]);
            d.set (1, y);
            c.add (d);
        }
        return c;
    }
    CubesVect mirrored () const
    {
        const int m = max_coord (0);
        CubesVect c;
        for (Dim d : cells)
        {
            d.set (0, m - d[0]);
            c.add (d);
        }
        return c;
    }
    bool maybe_canonical () const
    {
        const Box bb = bounding_box ();
        assert (bb.lo == Dim::all(0));
        for (int j = 1; j < bb.hi.size (); ++j)
            if (bb.hi[j] > bb.hi[j - 1])
                return false;
        return true;
    }
    bool is_canonical () const
    {
        return maybe_canonical () && *this == canonical ();
    }
    int multiplicity () const;
    const CubesVect& canonical () const;
    CubesVectIterator begin () const;
    CubesVectIterator end () const;
    CubesVectIterator cbegin () const;
    CubesVectIterator cend () const;
    friend CubesVectIterator;

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

struct Pad : std::vector<CubesVect>
{
    static Pad& get ();
    void add (const CubesVect &c)
    {
        for (auto it = begin (); ; ++it)
            if (int i; it == end () || (i = c.cmp (*it)) < 0)
            {
                std::vector<CubesVect>::insert (it, c);
                break;
            }
            else if (i == 0)
                break;
    }
};

using ScratchPads = std::vector<Pad>;
inline ScratchPads pads = ScratchPads (omp_get_max_threads ());

inline Pad& Pad::get ()
{
    return pads.at (omp_get_thread_num ());
}

const CubesVect& CubesVect::canonical () const
{
    CubesVect c (*this);
    Pad& pad = Pad::get ();
    pad.reserve (8);
    pad[0] = c;
    for (int i = 1; i < 8; ++i)
        if (c = i % 4 == 3 ? c.mirrored() : c.rotated(); c.maybe_canonical ())
            pad.add (c);
    return pad[0];
}

int CubesVect::multiplicity () const
{
    CubesVect c (*this);
    Pad& pad = Pad::get ();
    pad.reserve (8);
    pad[0] = c;
    for (int i = 1; i < 8; ++i)
        pad.add (c = i % 4 == 3 ? c.mirrored() : c.rotated());
    return pad.size ();
}

#undef  CUBES
#define CUBES CubesVect
#include "cubes-common.def"

#endif // CUBES_VECT_H

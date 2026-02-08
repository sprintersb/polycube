// -*- c++ -*-
#ifndef CUBES_VECT_H
#define CUBES_VECT_H

// Cubes is a vect since that is most memory friendly.

#include <string>
#include <vector>
// Own
#include "hash.h"
#include "dim.h"

#ifndef CUBES_H
#error use #include "cubes.h"
#endif

struct CubesIterator;

struct Cubes
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
    hash_t hash () const
    {
        hash_t h = CCITT32::crc32_init;
        for (Dim d : cells)
            h = CCITT32::crc32_add (h, d.ival());
        return h;
    }
    bool contains (Dim d) const
    {
        for (Dim c : cells)
        {
            const int i = c.cmp (d);
            if (i == 0)
                return true;
            if (i > 0)
                break;
        }
        return false;
    }
private:
    void add (Dim d)
    {
        for (auto it = cells.begin (); ; ++it)
        {
            int i;
            if (it == cells.end () || (i = d.cmp (*it)) < 0)
            {
                cells.insert (it, d);
                break;
            }
            if (i == 0)
                assert (i != 0 && "Assume we always increase cells");
        }
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
public:
    CubesIterator begin () const;
    CubesIterator end () const;
    CubesIterator cbegin () const;
    CubesIterator cend () const;
    friend CubesIterator;

#include "cubes-common.def"
};

struct CubesIterator
{
    Cubes::cell_iterator it;
    void operator ++ () { ++it; }
    bool operator == (const CubesIterator &r) const { return it == r.it; }
    bool operator != (const CubesIterator &r) const { return it != r.it; }
    Dim operator * () const { return *it; }
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

#endif // CUBES_VECT_H

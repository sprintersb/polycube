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
struct Pad;
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
    Pad congruents () const;
    bool is_canonical () const;
    void canonicalize ();
    int multiplicity () const;
    CubesVect canonical () const
    {
        CubesVect c (*this);
        c.canonicalize ();
        return c;
    }
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

Pad CubesVect::congruents () const
{
    using CV = CubesVect;
    using F = CV (*)(const CV&);
    CV c (*this);
    Pad pad;
    const char *S = "";
#if DIM == 2
    pad.reserve (8);
    static const F fun[] =
        {
            [](const CV &c) -> CV { return c.mirrored (0); },
            [](const CV &c) -> CV { return c.rotated (0, 1); },
        };
    S = "111" "0" "111";
#elif DIM == 3
    pad.reserve (48);
    static const F fun[] =
        {
            [](const CV &c) -> CV { return c.mirrored (0); },
            [](const CV &c) -> CV { return c.rotated (0, 1); },
            [](const CV &c) -> CV { return c.rotated (1, 2); },
            [](const CV &c) -> CV { return c.rotated (2, 0); },
        };
    S = ("12121 3 12121 3 12121 1 12121 0"
         "12121 3 12121 3 12121 1 12121");
#elif DIM == 4
    pad.reserve (2 * 192);
    static const F fun[] =
        {
            [](const CV &c) -> CV { return c.mirrored (0); },
            [](const CV &c) -> CV { return c.rotated (0, 1); }, // 1
            [](const CV &c) -> CV { return c.rotated (0, 2); }, // 2
            [](const CV &c) -> CV { return c.rotated (1, 2); }, // 3
            [](const CV &c) -> CV { return c.rotated (2, 1); }, // 4
            [](const CV& c) -> CV { return c.rotated (3, 0); }, // 5
            [](const CV& c) -> CV { return c.rotated (3, 1); }, // 6
        };
#define M "11121113111311141112111"
    S = (M "5" M "5" M "6" M "5" M "6" M "5" M "5" M "0"
         M "5" M "5" M "6" M "5" M "6" M "5" M "5" M);
#undef M
#else
    static const F fun[] = {};
    assert (0 && "todo: canonicalize in DIM");
#endif
    pad.resize (1, c);
    for (auto s = S; *s; ++s)
        if (*s != ' ')
            pad.add (c = fun[*s - '0'] (c));
    return pad;
}

void CubesVect::canonicalize ()
{
    cells = std::move (congruents()[0].cells);
}

bool CubesVect::is_canonical () const
{
    return *this == congruents()[0];
}

int CubesVect::multiplicity () const
{
    return (int) congruents().size ();
}

#undef  CUBES
#define CUBES CubesVect
#include "cubes-common.def"

#endif // CUBES_VECT_H

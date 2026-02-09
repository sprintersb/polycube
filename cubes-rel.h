// -*- c++ -*-
#ifndef CUBES_REL_H
#define CUBES_REL_H

// CubesRel refers to smaller cubes for memory footprint.

#include <string>
#include <set>

#include "hash.h"
#include "dim.h"

#ifndef CUBES_H
#error use #include "cubes.h"
#endif

struct CubesRel
{
    struct CIterator;

    const CubesRel *dad = nullptr;
    Dim cube;

    CubesRel () {}
    CubesRel (const CubesRel *dad, Dim d) : dad(dad), cube(d) {}

    int size () const
    {
        int sz = 0;
        for (Dim d : *this)
            ++sz, (void) d;
        return sz;
    }
    Dim min_cube () const
    {
        Dim d (cube);
        for (auto p = dad; p; p = p->dad)
            d.min (p->cube);
        return d;
    }
    // As is seems, for the sizes of interest, the execution times are:
    // CubesVect(array) < CubesVect(vector) < set<Dim> < unordered_set<Dim>.
    using SortedCubes = CubesVect;
    SortedCubes sorted () const
    {
        const Dim m = min_cube ();
        SortedCubes set;
        for (Dim d : *this)
            set.add (d - m);
        return set;
    }
    // Symmetric in cubes and shift-invariant but quite expensive.
    bool operator == (const CubesRel &c) const
    {
        return sorted () == c.sorted ();
    }
    bool is_canonical () const
    {
        return sorted().is_canonical ();
    }
    int multiplicity () const
    {
        return sorted().multiplicity ();
    }
    // Shift-invariant and symmetric.
    Hasher::type hash () const
    {
        const Dim m = min_cube ();
        Hasher::type h = 0;
        for (Dim d : *this)
            h += Hasher() ((d - m).ival());
        return h;
    }
    bool contains (Dim d) const
    {
        for (Dim c : *this)
            if (c == d)
                return true;
        return false;
    }

    struct CIterator
    {
        const CubesRel *ptr;
        void operator ++ () { ptr = ptr->dad; }
        bool operator == (const CIterator &r) const { return ptr == r.ptr; }
        bool operator != (const CIterator &r) const { return ptr != r.ptr; }
        Dim operator * () const { return ptr->cube; }
    };
    using const_iterator = CIterator;
    CIterator begin () const { return CIterator { this }; }
    CIterator end ()   const { return CIterator { nullptr }; }
    CIterator cbegin () const { return CIterator { this }; }
    CIterator cend ()   const { return CIterator { nullptr }; }

    // Common tail.
public:
    Box bounding_box () const;
    std::string ascii (char c = '*') const;
    friend std::ostream& operator << (std::ostream&, const CubesRel&);
};

#undef  CUBES
#define CUBES CubesRel
#include "cubes-common.def"

#endif // CUBES_REL_H

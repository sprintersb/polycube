// -*- c++ -*-
#ifndef CUBES_REL_H
#define CUBES_REL_H

// Cubes refers to smaller cubes for memory footprint.

#include <string>
#include <set>

#include "hash.h"
#include "dim.h"

#ifndef CUBES_H
#error use #include "cubes.h"
#endif

struct Cubes
{
    struct CIterator;

    const Cubes *dad = nullptr;
    Dim cube;

    Cubes () {}
    Cubes (const Cubes *dad, Dim d) : dad(dad), cube(d) {}

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
    using SortedCubes = std::set<Dim>;
    SortedCubes sorted () const
    {
        const Dim m = min_cube ();
        SortedCubes set;
        for (Dim d : *this)
            set.emplace (d - m);
        return set;
    }
    // Symmetric in cubes and shift-invariant but quite expensive.
    bool operator == (const Cubes &c) const
    {
        return sorted () == c.sorted ();
    }
    // Shift-invariant and symmetric.
    hash_t hash () const
    {
        const Dim m = min_cube ();
        hash_t h = 0;
        for (Dim d : *this)
            h += CCITT32::crc32_value ((d - m).ival());
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
        const Cubes *ptr;
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

#include "cubes-common.def"
};

using CubesIterator = Cubes::CIterator;

#endif // CUBES_REL_H

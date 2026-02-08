#ifndef CUBES_H
#define CUBES_H

#ifdef CUBES_REL
#undef CUBES_ARRAY
#endif

#if defined CUBES_ARRAY && !defined CELLS
#error CUBES_ARRAY CELLS=?
#endif

#include <iostream>
#include "dim.h"

class Cubes;

#ifdef CUBES_REL
#include "cubes-rel.h"
#else
#include "cubes-vect.h"
#endif

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

inline std::ostream& operator << (std::ostream &ost, const Cubes &c)
{
    ost << "{#" << c.size ();
    for (Dim d : c)
        ost << " " << d;
    return ost << " }";
}

#include "cubes-border.h"

#endif // CUBES_H

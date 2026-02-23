// -*- c++ -*-
#ifndef CUBES_BORDER_H
#define CUBES_BORDER_H

#include <iostream>
#include <string>
#include <list>
#include <vector>
// Own
#include "cubes.h"

struct Polygon : public std::list<Dim>
{
    std::string svg (bool rel = 1) const;
    void push (Dim d, bool tidy);
    friend std::ostream& operator << (std::ostream&, const Polygon&);
};

struct Polygons : public std::vector<Polygon>
{
    std::string svg (bool rel = 1) const
    {
        std::string s;
        for (const Polygon &p : *this)
            s += std::string (" "  + !s.size()) + p.svg (rel);
        return s;
    }
    friend std::ostream& operator << (std::ostream&, const Polygons&);
};

extern Polygons cubes_border (const Cubes&);

#endif // CUBES_BORDER_H

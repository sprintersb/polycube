// -*- c++ -*-
#ifndef CUBES_BORDER_H
#define CUBES_BORDER_H

#include <utility> // std::move, std::swap
#include <iostream>
#include <sstream>
#include <string>

#include <list>
#include <unordered_set>
#include <unordered_map>
// C
#include <cassert>
// Own
#include "hash.h"
#include "dim.h"

// Lines are only used in dimension 2, so they occupy only 4 bytes.
// Hence pass them as objects and not as const Line&.
struct Line
{
    Dim a, b;
    bool operator == (Line l) const
    {
        return a == l.a && b == l.b;
    }
    struct Hash
    {
        hash_t operator () (Line l) const
        {
            return l.a.hash() + 13 * l.b.hash();
        }
    };
    struct SymmetricHash
    {
        hash_t operator () (Line l) const
        {
            return l.a.hash() + l.b.hash();
        }
    };
    struct SymmetricEQ
    {
        bool operator () (Line g, Line h) const
        {
            return g == h || (g.a == h.b && g.b == h.a);
        }
    };
    friend std::ostream& operator << (std::ostream &ost, Line l);
};

struct Polygon : public std::list<Dim>
{
    std::string svg (bool rel = 1) const
    {
        std::stringstream ss;
        Dim from{};
        const Dim O;
        bool start = true;
        for (const Dim d : *this)
        {
            if (start)
            {
                ss << "M" << d[0] << "," << d[1];
                start = false;
            }
            else if ((d - from)[0] == 0)
                ss << "Vv"[rel] << (d - (rel ? from : O))[1];
            else if ((d - from)[1] == 0)
                ss << "Hh"[rel] << (d - (rel ? from : O))[0];
            else
                assert (0);
            from = d;
        }
        ss << "Zz"[rel];
        return ss.str ();
    }
    void push (Dim d, bool tidy)
    {
        if (auto p0 = rbegin (); tidy && p0 != rend ())
            if (auto p1 = p0; ++p1 != rend ())
            {
                Dim d1 = *p0 - *p1;
                Dim d2 = d - *p0;
                if (d1[0] == 0 && d2[0] == 0)
                {
                    p0->v[1] = d[1];
                    return;
                }
                else if (d1[1] == 0 && d2[1] == 0)
                {
                    p0->v[0] = d[0];
                    return;
                }
            }
        push_back (d);
    }
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

struct BorderFinder
{
    struct LinePool
    {
        using LineBucket = std::unordered_set
            <Line, Line::SymmetricHash, Line::SymmetricEQ>;
        std::unordered_map<Dim, LineBucket, Dim::Hash> bucks;
        void operator += (Line l)
        {
            add (l.a, l);
            add (l.b, l);
        }
        int operator -= (Line l)
        {
            bool killa = sub(l.a, l);
            bool killb = sub(l.b, l);
            return killa + killb;
        }
        void add (Dim d, Line l)
        {
            if (auto &&b = bucks.find (d); b == bucks.end ())
                bucks[d] = LineBucket { l };
            else
                b->second.emplace (l);
        }
        bool sub (Dim d, Line l)
        {
            if (auto &&b = bucks.find (d); b == bucks.end ())
                return false;
            else
            {
                const bool erased = b->second.erase (l);
                if (b->second.empty ())
                    bucks.erase (d);
                return erased;
            }
        }
        const LineBucket* get (Dim d) const
        {
            auto &&b = bucks.find (d);
            return b == bucks.end () ? nullptr : & b->second;
        }
        bool contains (Line line) const
        {
            const LineBucket *buck = get (line.a);
            return buck && buck->find (line) != buck->end ();
        }
        friend std::ostream& operator << (std::ostream&, const LinePool&);
    }; // LinePool

    const Cubes &cs;
    LinePool pool;
    Polygons polygons;

    BorderFinder (const Cubes &cs)
        : cs(cs)
    {
        assert (DIM == 2);
        // Collect all border line segments in pool.
        for (Dim p : cs)
            for (Dim d : p)
                if (! cs.contains (p + d))
                {
                    // P+d is in corona.  Determine unoriented line
                    // segment between P and P+d.
                    Dim off = Dim { d[0] + d[1] > 0, d[1] > d[0] };
                    pool += Line { p + off, p + off + d.rot(1) };
                }
    }
    enum Orient { Left = 1, Right = -1 };
    Line get_start (bool outer) const
    {
        if (outer)
        {
            for (Dim p : cs)
                if (p[1] == 0)
                    return { p, p + Dim{1,0} };
        }
        else
            for (Dim p : cs)
                for (Dim::value_t dy = 0; dy <= 1; ++dy)
                {
                    const Dim a = p + Dim { 0, dy };
                    const Dim b = a + Dim { 1, 0 };
                    Line line { dy ? b : a, dy ? a : b };
                    if (pool.contains (line))
                        return line;
                }
        assert (0);
    }
    Polygon get_polygon (bool outer, bool tidy)
    {
        Polygon polygon;
        const Line start = get_start (outer);
        Orient orient = start.b[0] > start.a[0] ? Left : Right;
        Dim p2{};
        for (Line line = start; ; line = Line { line.b, p2 })
        {
            polygon.push (line.a, tidy);
            assert (2 == (pool -= line));
            if (line.b == start.a)
                break;

            // Iterate over all pool lines that start / end at line.b.
            const LinePool::LineBucket *buck = pool.get (line.b);
            assert (buck && ! buck->empty ());
            const Dim dir = line.b - line.a;
            int sin_angle = -100;
            for (Line l : *buck)
            {
                if (l.a != line.b)
                    std::swap (l.a, l.b);
                assert (l.a == line.b);
                const int sina = (dir % (l.b - l.a)) * orient;
                // > makes the outer border greedy, < makes it reluctant.
                if (sina > sin_angle)
                    p2 = l.b, sin_angle = sina;
            }
            assert (sin_angle >= -1 && sin_angle <= 1);
        }
        polygon.push (start.a, tidy);
        return polygon;
    }
    // Return a vector of Polygons that represent the border of the
    // Cubes.  The first polygon is the outer border and maybe more.
    // The rest represent (parts of) the inner border and is optional.
    Polygons border (bool tidy = 1)
    {
        for (bool outer = 1; ! pool.bucks.empty (); outer = 0)
            polygons.push_back (get_polygon (outer, tidy));
        assert (polygons.size () >= 1);
        return polygons;
    }
};

inline std::ostream& operator << (std::ostream &ost, Line l)
{
    return ost << l.a << "--" << l.b;
}

inline std::ostream& operator << (std::ostream &ost,
                                  const BorderFinder::LinePool &lp)
{
    for (const auto& b : lp.bucks)
    {
        ost << b.first << ":";
        for (auto l : b.second)
            ost << "  " << l;
        ost << "\n";
    }
    return ost;
}

inline std::ostream& operator << (std::ostream &ost, const Polygon &pg)
{
    ost << "Polygon: ";
    for (Dim p : pg)
        ost << " " << p;
    return ost << "\n";
}

#endif // CUBES_BORDER_H

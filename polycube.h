// -*- c++ -*-
#include <vector>
#include <list>
#include <set>
#include <map>
#include <utility> // std::move
#include <iostream>
#include <mutex>

#include <cstdint>
#include <cassert>
#include <omp.h>

class Dim;
class Cells;
class PolyCube;
inline std::ostream& operator << (std::ostream&, const Dim&);
inline std::ostream& operator << (std::ostream&, const Cells&);
inline std::ostream& operator << (std::ostream&, const PolyCube&);

struct Dim
{
    using value_type = int8_t;
    std::vector<value_type> x;
    Dim (const std::vector<value_type> &x) : x(x) {}
    Dim (std::vector<value_type> &&x) : x(std::move (x)) {}
    static Dim dim (int n)
    {
        return std::vector<value_type> (n, 0);
    }
    int size () const
    {
        return (int) x.size ();
    }
    int cmp (const Dim &d) const
    {
        if (size () != d.size ())
            return size () - d.size ();
        for (int i = 0; i < size (); ++i)
            if (x[i] != d.x[i])
                return x[i] - d.x[i];
        return 0;
    }
};

struct Cells
{
    std::list<Dim> cells;

    int size () const
    {
        return (int) cells.size ();
    }
    void add (const Dim &d)
    {
        auto end = cells.end ();
        for (auto p = cells.begin(); ; ++p)
            if (p == end)
            {
                cells.emplace (p, d);
                break;
            }
            else
            {
                const int i = d.cmp (*p);
                if (i > 0)
                    continue;
                if (i < 0)
                    cells.emplace (p, d);
                break;
            }
    }
    void erase (const Dim &d)
    {
        auto end = cells.end ();
        for (auto p = cells.begin(); p != end; ++p)
        {
            const int i = p->cmp (d);
            if (i == 0)
                cells.erase (p);
            if (i >= 0)
                break;
        }
    }
    int cmp (const Cells &c) const
    {
        auto p2 = c.cells.begin ();
        auto e2 = c.cells.end ();
        for (const auto &c : cells)
        {
            if (p2 == e2)
                return 1;
            const int i = c.cmp (*p2);
            if (i)
                return i;
            ++p2;
        }
        return p2 == e2 ? 0 : -1;
    }
    bool contains (const Dim &d) const
    {
        for (const auto &c : cells)
        {
            const int i = c.cmp (d);
            if (i == 0)
                return true;
            if (i > 0)
                break;
        }
        return false;
    }
    void shift (int i, int off)
    {
        for (auto &c : cells)
        {
            assert (i < c.size ());
            c.x[i] += off;
        }
    }
};

struct PolyCube
{
    using Set = std::set<PolyCube>;
    using List = std::list<PolyCube>;
    struct MuxSet
    {
        // std::mutex is not movable, so we have to use the
        // swap trick to get a vector of given size, see
        // https://stackoverflow.com/a/24170141/1556746
        std::mutex mux;
        Set set;
    };
    using Vector = std::vector<MuxSet>; // Indexed by corona size.
    Cells cubes;
    Cells corona;
    void add (const Dim &d)
    {
        cubes.add (d);
        corona.erase (d);
        // Repair corona.
        for (int i = 0; i < d.size (); ++i)
            for (int dir = 1; dir >= -1; dir -= 2)
            {
                Dim d2 (d);
                d2.x[i] += dir;
                if (! cubes.contains (d2))
                    corona.add (d2);
            }
        // Normalize cubes.
        for (int i = 0; i < d.size (); ++i)
            if (d.x[i] < 0)
            {
                cubes.shift (i, -d.x[i]);
                corona.shift (i, -d.x[i]);
            }
    }
    bool operator < (const PolyCube &c) const
    {
        return cubes.cmp (c.cubes) < 0;
    }
    // Way 0
    void add_sprouts (Set &set) const
    {
        for (const Dim &d : corona.cells)
        {
            PolyCube pc (*this);
            pc.add (d);
            set.emplace (std::move (pc));
        }
    }

    // Way 1
    void add_sprouts_omp (Set &set) const
    {
        std::vector<const Dim*> dim;
        dim.resize (corona.size ());
        int i = 0;
        for (const Dim &d : corona.cells)
            dim[i++] = &d;

#pragma omp parallel for
        for (int i = 0; i < corona.size (); ++i)
        {
            PolyCube pc (*this);
            pc.add (* dim[i]);
#pragma omp critical
            set.emplace (std::move (pc));
        }
    }

    // Way 2
    void add_sprouts (List &list) const
    {
        for (const Dim &d : corona.cells)
        {
            PolyCube pc (*this);
            pc.add (d);
            list.emplace_back (std::move (pc));
        }
    }

    // Way 2
    static void add_sprouts (Set &set2, const Set &set)
    {
        std::vector<const PolyCube*> pc;
        pc.resize (set.size ());
        int j = 0;
        for (const auto &p : set)
            pc[j++] = &p;

#pragma omp parallel for
        for (size_t j = 0; j < pc.size (); ++j)
        {
            List list;
            pc[j]->add_sprouts (list);
#pragma omp critical
            for (auto &p : list)
                set2.emplace (std::move (p));
        }
    }

    // Way 3
    void add_sprouts (Vector &vms) const
    {
        for (const Dim &d : corona.cells)
        {
            PolyCube pc (*this);
            pc.add (d);
            auto &slot = vms[pc.corona.size ()];
            slot.mux.lock ();
            slot.set.emplace (std::move (pc));
            slot.mux.unlock ();
        }
    }

    // Way 3
    static void add_sprouts (int dim, int n, Vector &vset2, const Vector &vset)
    {
        assert (n >= 2 && dim >= 1);
        const int max_corona_size = 2 * (dim - 1) * n + 2;

        Vector v (1 + max_corona_size);
        vset2.swap (v); // Since resize() doesn't like std::mutex.

        size_t n_polycubes = 0;
        for (const auto &ms : vset)
            n_polycubes += ms.set.size ();
        std::vector<const PolyCube*> vpc (n_polycubes);
        int j = 0;
        for (const auto &ms : vset)
            for (const auto &pc : ms.set)
                vpc[j++] = &pc;
#pragma omp parallel for
        for (size_t j = 0; j < vpc.size (); ++j)
            vpc[j]->add_sprouts (vset2);
    }

    // Univariate polynomial over Z in sparse representation.
    using Poly = std::map<int,int>;
    static Poly get_poly (const Set &set)
    {
        Poly poly;
        for (const auto &pc : set)
        {
            const int coro = pc.corona.size ();
            const auto &monome = poly.find (coro);
            if (monome == poly.end ())
                poly[coro] = 1;
            else
                monome->second += 1;
        }
        return poly;
    }

    static Poly get_poly (const Vector &vms)
    {
        Poly poly;
        for (size_t j = 0; j < vms.size (); ++j)
            if (vms[j].set.size ())
                poly[j] = vms[j].set.size ();
        return poly;
    }

    static void print_poly (int n, const Poly &poly, const char *var = "q")
    {
        printf ("c");
        printf (n > 9 ? "_{%d}(p) = p" : "_%d(p) = p", n);
        if (n != 1) printf (n > 9 ? "^{%d}" : "^%d", n);
        printf (" \\cdot ");
        printf (poly.size() > 1 ? "(" : "");
        bool start = true;
        for (auto m : poly)
        {
            if (!start)
                printf (" + ");
            if (m.second != 1)
                printf ("%d", m.second);
            printf (m.first < 10 ? "%s^%d" : "%s^{%d}", var, m.first);
            start = false;
        }
        printf ("%s\n", poly.size() > 1 ? ")" : "");

        start = true;
        printf ("/*%d*/ { ", n);
        for (auto m : poly)
        {
            if (!start)
                printf (", ");
            printf ("%d,%d", m.second, m.first);
            start = false;
        }
        printf (" }\n");
    }
};

inline std::ostream& operator << (std::ostream &ost, const PolyCube &pc)
{
    ost << "cubes: " << pc.cubes << "\n";
    ost << "coron: " << pc.corona << "\n";
    return ost;
}

inline std::ostream& operator << (std::ostream &ost, const Dim &d)
{
    ost << "<";
    const char *comma = "";
    for (auto x : d.x)
    {
        ost << comma << (int) x;
        comma = ",";
    }
    return ost << ">";
}

inline std::ostream& operator << (std::ostream &ost, const Cells &c)
{
    ost << "{#" << c.size ();
    for (auto c : c.cells)
        ost << " " << c;
    return ost << " }";
}

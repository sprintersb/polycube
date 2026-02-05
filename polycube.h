// -*- c++ -*-
#include <utility> // std::move
#include <iostream>
#include <mutex>
#include <atomic>
#include <string>
#include <sstream>
// Containers
#include <array>
#include <list>
#include <map>
#include <vector>
#include <set>
#include <unordered_set>
#include <unordered_map>
#include <iterator>
// C
#include <cstdint>
#include <cinttypes>
#include <cassert>
#include <climits>
// Other
#include <omp.h>
// Own
#include "polycube-count.h"

class Dim;
class Cubes;
class Corona;
class PolyCube;
inline std::ostream& operator << (std::ostream&, Dim);
inline std::ostream& operator << (std::ostream&, const Cubes&);
inline std::ostream& operator << (std::ostream&, const Corona&);
inline std::ostream& operator << (std::ostream&, const PolyCube&);

struct DimIterator;
using hash_t = uint32_t;

struct Dim
{
    using value_t = int8_t;
#if DIM == 2
    using vector_t = value_t __attribute__((vector_size(2)));
    using int_t = uint16_t;
#elif DIM > 2 && DIM <= 4
    using vector_t = value_t __attribute__((vector_size(4)));
    using int_t = uint32_t;
#elif DIM > 4 && DIM <= 8
    using vector_t = value_t __attribute__((vector_size(8)));
    using int_t = uint64_t;
#else
    using vector_t = value_t __attribute__((vector_size(16)));
    using int_t = unsigned __int128;
#endif

    vector_t v = (vector_t) (int_t) 0;
    static inline constexpr vector_t all0 = (vector_t) (int_t) 0;
    Dim () : v(Dim::all0) {}
    Dim (std::initializer_list<value_t> &&il)
    {
        assert (il.size () == 0 || il.size () == DIM);
        v = Dim::all0;
        int j = 0;
        for (auto i : il)
            v[j++] = i;
    }
    Dim (const vector_t &v) : v(v) {}

    static Dim all (int w)
    {
        Dim d{};
        for (int i = 0; i < d.size (); ++i)
            d.v[i] = w;
        return d;
    }
    int size () const
    {
        return DIM;
    }
    int operator [] (int i) const
    {
        return v[i];
    }
    // Shift-invariant comparison, so (int_t) <=> (int_t) won't work.
    // Benefit is that Cubes.cells don't change their order when shifted.
    int cmp (Dim d) const
    {
        for (int i = 0; i < size (); ++i)
            if (v[i] != d.v[i])
                return v[i] - d.v[i];
        return 0;
    }
    bool operator == (Dim d) const { return (int_t) v == (int_t) d.v; }
    bool operator != (Dim d) const { return (int_t) v != (int_t) d.v; }
    bool operator <= (Dim d) const { return cmp (d) <= 0; }
    bool operator >= (Dim d) const { return cmp (d) >= 0; }
    bool operator <  (Dim d) const { return cmp (d) <  0; }
    bool operator >  (Dim d) const { return cmp (d) >  0; }
    void operator += (Dim d) { v += d.v; }
    void operator -= (Dim d) { v -= d.v; }
    void operator *= (int i) { v *= (value_t) i; }
    Dim operator + (Dim d) const { return Dim (v + d.v); }
    Dim operator - (Dim d) const { return Dim (v - d.v); }
    Dim operator * (int i) const { return Dim (v * (value_t) i); }
    int operator % (Dim d) const { return v[0] * d[1] - v[1] * d[0]; }
    Dim rot (int i /* Left in units of 90 deg */) const
    {
        Dim d (*this);
        i = (4 + (i % 4)) % 4;
        while (i-- > 0)
            d = Dim { (value_t) -d.v[1], d.v[0] };
        return d;
    }
    void min (Dim d)
    {
        for (int i = 0; i < size (); ++i)
            v[i] = std::min (v[i], d.v[i]);
    }
    void max (Dim d)
    {
        for (int i = 0; i < size (); ++i)
            v[i] = std::max (v[i], d.v[i]);
    }
    hash_t hash () const
    {
        if (sizeof (int_t) <= sizeof (hash_t))
            return (hash_t) (int_t) v;
        else
        {
            hash_t h = 0;
            for (int i = 0; i < size (); ++i)
                h = 13 * h + (hash_t) v[i];
            return h;
        }
    }
    struct Hash
    {
        hash_t operator () (Dim d) const
        {
            return d.hash ();
        }
    };

    DimIterator begin () const;
    DimIterator end () const;
};

struct Box
{
    Dim lo, hi;
    bool contains (Dim d)
    {
        for (int i = 0; i < d.size (); ++i)
            if (d.v[i] < lo.v[i] || d.v[i] > hi.v[i])
                return false;
        return true;
    }
    Box grow (int g) const
    {
        return Box { lo - Dim::all (g),  hi + Dim::all (g) };
    }
};


class DimIterator
{
    using Corona0 = std::array<Dim, 2 * DIM>;
    static inline Corona0 get_corona0 ()
    {
        Corona0 c;
        for (int i = 0; i < 2 * DIM; ++i)
            c[i].v[i / 2] = i % 2 == 0 ? 1 : -1;
        return c;
    }
    static inline const Corona0 corona0 = DimIterator::get_corona0 ();
    int pos;

    DimIterator (int pos) : pos(pos) {}
    friend Dim;
public:
    bool operator != (DimIterator di) const
    {
        return pos != di.pos;
    }
    void operator ++ ()
    {
        ++ pos;
    }
    Dim operator * () const
    {
        return DimIterator::corona0[pos];
    }
};

inline DimIterator Dim::begin () const { return DimIterator (0); }
inline DimIterator Dim::end ()   const { return DimIterator (2 * DIM); }


// Is a vect since that is most memory friendly.
struct Cubes
{
    std::vector<Dim> cells;

    int size () const
    {
        return cells.size ();
    }
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
    }
    Cubes operator * (int i) const
    {
        Cubes c { *this };
        for (Dim &d : c.cells)
            d *= i;
        return c;
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
    hash_t hash () const
    {
        hash_t h = 0;
        for (Dim d : cells)
            h = h * 13 + d.hash ();
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
    void shift (int i, int off)
    {
        for (auto &c : cells)
        {
            assert (i < c.size ());
            c.v[i] += off;
        }
    }

    Box bounding_box () const
    {
        Box box { cells[0], cells[0] };
        for (int i = 1; i < size (); ++i)
        {
            box.lo.min (cells[i]);
            box.hi.max (cells[i]);
        }
        return box;
    }

    //////////////////////////////////////////////////////////////////////////
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
                    ss << "M" << d[0] << " " << d[1];
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
    using Polygons = std::vector<Polygon>;

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
            for (Dim p : cs.cells)
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
                for (Dim p : cs.cells)
                    if (p[1] == 0)
                        return { p, p + Dim{1,0} };
            }
            else
                for (Dim p : cs.cells)
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
    //////////////////////////////////////////////////////////////////////////

    std::string ascii (char c = '*') const
    {
#if DIM == 2
        auto bbox = bounding_box ();
        std::string str;
        for (Dim::value_t y = bbox.hi[1]; y >= bbox.lo[1]; --y)
        {
            for (Dim::value_t x = bbox.lo[0]; x <= bbox.hi[0]; ++x)
                str += contains (Dim { x, y }) ? c : ' ';
            str += "\n";
        }
        return str;
#else
        (void) c;
        return "Cubes.ascii(DIM != 2)";
#endif
    }
};


struct Corona
{
    std::unordered_set<Dim, Dim::Hash> cells;

    int size () const
    {
        return cells.size ();
    }
    void add (Dim d)
    {
        cells.emplace (d);
    }
    void erase (Dim d)
    {
        cells.erase (d);
    }
    bool contains (Dim d) const
    {
        return cells.find (d) != cells.end ();
    }
    Corona operator * (int i) const
    {
        Corona c;
        for (Dim d : cells)
            c.add (d * i);
        return c;
    }
};

#define BACK "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b"


struct PolyCube
{
    struct Hash
    {
        hash_t operator () (const PolyCube &pc) const
        {
            return pc.hash ();
        }
    };

    using Set = std::unordered_set<PolyCube, Hash>;

#pragma omp declare reduction(merge : Set : merge (omp_out, omp_in)) \
    initializer (omp_priv = omp_orig)
    static void merge (Set &sout, Set &sin)
    {
        sout.merge (sin);
    }
    static void merge (Set &s, PolyCube &&pc)
    {
        s.emplace (std::move (pc));
    }

    using List = std::list<PolyCube>;
    struct MuxSet
    {
        // std::mutex is not movable, so we have to use the
        // swap trick to get a vector of given size, see
        // https://stackoverflow.com/a/24170141/1556746
        std::mutex mux;
        Set set;
    };
    using Vector = std::vector<MuxSet>;
    Cubes m_cubes;
    hash_t m_hash = 0;

    bool contains (Dim d) const
    {
        return m_cubes.contains (d);
    }
    Corona corona () const
    {
        Corona cora;
        for (Dim d : m_cubes.cells)
            for (Dim delta : d)
                if (! m_cubes.contains (d + delta))
                    cora.add (d + delta);
        return cora;
    }
    Box bounding_box () const
    {
        return m_cubes.bounding_box ();
    }
    PolyCube operator * (int i) const
    {
        Cubes c = m_cubes * i;
        return PolyCube { c, c.hash () };
    }
    bool has_large_corona (int max_corona) const
    {
        const Corona &cora = corona ();
        if (cora.size () > max_corona)
            return true;
        // Try some simple convexity tests.  Use an unordered_set for faster
        // accesses.  Additional size doesn't matter here.
        Corona cubs;
        for (const auto &c : m_cubes.cells)
            cubs.add (c);
        for (Dim d : cora.cells)
        {
            if (cubs.contains (d + Dim{1,0}) && cubs.contains (d - Dim{1,0}))
                return true;
            if (cubs.contains (d + Dim{0,1}) && cubs.contains (d - Dim{0,1}))
                return true;
        }
        return false; // Only weakly false, i.e. not necessarily convex.
    }
    void add (Dim d)
    {
        m_cubes.add (d);
        // Normalize m_cubes.
        for (int i = 0; i < d.size (); ++i)
            if (d.v[i] < 0)
                m_cubes.shift (i, -d.v[i]);
        m_hash = m_cubes.hash ();
    }
    bool operator == (const PolyCube &c) const
    {
        return m_cubes.cmp (c.m_cubes) == 0;
    }
    bool operator < (const PolyCube &c) const
    {
        return m_cubes.cmp (c.m_cubes) < 0;
    }
    hash_t hash () const
    {
        return m_hash;
    }
    static Set find_min_corona (const Vector &vms)
    {
        Set set;
        std::atomic<int> corona_size = INT_MAX;
#pragma omp parallel for schedule(dynamic) // reduction(min: corona_size)
        for (size_t j = 0; j < vms.size (); ++j)
            for (const auto &pc : vms[j].set)
            {
                int csize = pc.corona ().size();
#pragma omp critical
                if (csize <= corona_size)
                {
                    if (csize < corona_size)
                        set.clear ();
                    corona_size = csize;
                    set.insert (pc);
                }
            }
        return set;
    }

    // Way 0
    void add_sprouts (Set &set) const
    {
        for (Dim d : corona().cells)
        {
            PolyCube pc (*this);
            pc.add (d);
            set.emplace (std::move (pc));
        }
    }

    // Way 4
    int add_sprouts_way4 (Vector &vms, int max_corona) const
    {
        int new_count = 0;
        for (Dim d : corona().cells)
        {
            PolyCube pc (*this);
            pc.add (d);
            if (max_corona > 0 && pc.has_large_corona (max_corona))
                continue;
            MuxSet &slot = vms[pc.hash () % vms.size ()];

            slot.mux.lock ();
            const auto n = slot.set.size ();
            slot.set.emplace (std::move (pc));
            new_count += n != slot.set.size ();
            slot.mux.unlock ();
        }
        return new_count;
    }

    // Way 4
    static void add_sprouts_way4 (int n_cells, int n_slots,
                                  Vector &vset2, const Vector &vset,
                                  int max_corona = -1)
    {
        Vector v (n_slots);
        vset2.swap (v); // Since resize() doesn't like std::mutex

        // Only for printing stat.
        const int64_t n_cubes = cube_count (DIM, n_cells);
        std::atomic<int64_t> pc_count = 0;
        int64_t pc_show = pc_count;

#pragma omp parallel for schedule(dynamic)
        for (size_t j = 0; j < vset.size (); ++j)
        {
            for (const auto &pc : vset[j].set)
                pc_count += pc.add_sprouts_way4 (vset2, max_corona);

            // Print stat.
            if (omp_get_thread_num () == 0)
            {
                const int64_t pcs = pc_count;
                if (pc_show == 0 || pcs - pc_show > 100000)
                {
                    pc_show = pcs;
                    printf ("%s%" PRIi64 " Cubs", BACK, pcs);
                    if (n_cubes > 1)
                        printf (" = %.1f%%", 100.0 * pcs / n_cubes);
                    fflush (stdout);
                }
            } // master
        } // parallel for

        printf ("%s                                    %s%s", BACK, BACK, BACK);
        fflush (stdout);
    }

    // Univariate polynomial over Z in sparse representation.
    struct Poly
    {
        enum { POLY_TEX, POLY_LIST };
        std::map<int,int> a_;

        Poly () {}

        // Way 0: 100% sequential.
        Poly (const Set &set)
        {
            for (const auto &pc : set)
            {
                const int coro = pc.corona().size ();
                const auto &monome = a_.find (coro);
                if (monome == a_.end ())
                    a_[coro] = 1;
                else
                    monome->second += 1;
            }
        }

        // Way 4.
        Poly (const Vector &vms)
        {
            Poly::init (*this, vms);
        }

        // Way 4.
        Poly& operator += (const Poly &q)
        {
            for (const auto &q_mono : q.a_)
            {
                const int expo = q_mono.first;
                const int coef = q_mono.second;
                const auto &p_mono = a_.find (expo);
                if (p_mono == a_.end ())
                    a_[expo] = coef;
                else
                    p_mono->second += coef;
            }
            return *this;
        }

#pragma omp declare                                                     \
    reduction (+ : Poly : omp_out += omp_in)                            \
    initializer (omp_priv = omp_orig)

        static void init (Poly &poly, const Vector &vms)
        {
            // Way 4.
#pragma omp parallel for schedule(dynamic,20) reduction(+: poly)
            for (size_t j = 0; j < vms.size (); ++j)
                poly += Poly (vms[j].set);
        }

        void print (int n, int style, const char *var = "q") const
        {
            if (style == POLY_TEX)
            {
                printf ("c");
                printf (n > 9 ? "_{%d}(p) = p" : "_%d(p) = p", n);
                if (n != 1) printf (n > 9 ? "^{%d}" : "^%d", n);
                printf (" \\cdot ");
                printf (a_.size() > 1 ? "(" : "");
                bool start = true;
                for (auto m : a_)
                {
                    if (!start)
                        printf (" + ");
                    if (m.second != 1)
                        printf ("%d", m.second);
                    printf (m.first < 10 ? "%s^%d" : "%s^{%d}", var, m.first);
                    start = false;
                }
                printf ("%s\n", a_.size() > 1 ? ")" : "");
            }
            else if (style == POLY_LIST)
            {
                bool start = true;
                printf ("/*%d*/ { ", n);
                for (auto m : a_)
                {
                    if (!start)
                        printf (", ");
                    printf ("%d,%d", m.second, m.first);
                    start = false;
                }
                printf (" }\n");
            }
        }
    }; // Poly
}; // PolyCube

inline std::ostream& operator << (std::ostream &ost, const PolyCube &pc)
{
    ost << "cubes: " << pc.m_cubes << "\n";
    ost << "coron: " << pc.corona() << "\n";
    return ost;
}

inline std::ostream& operator << (std::ostream &ost, Dim d)
{
    char comma = '<';
    for (int i = 0; i < d.size (); ++i)
    {
        ost << comma << (int) d[i];
        comma = ',';
    }
    return ost << '>';
}

inline std::ostream& operator << (std::ostream &ost, const Corona &c)
{
    ost << "{#" << c.size ();
    for (auto c : c.cells)
        ost << " " << c;
    return ost << " }";
}

inline std::ostream& operator << (std::ostream &ost, const Cubes &c)
{
    ost << "{#" << c.size ();
    for (auto c : c.cells)
        ost << " " << c;
    return ost << " }";
}

inline std::ostream& operator << (std::ostream &ost, Cubes::Line l)
{
    return ost << l.a << "--" << l.b;
}

inline std::ostream& operator << (std::ostream &ost,
                                  const Cubes::BorderFinder::LinePool &lp)
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

inline std::ostream& operator << (std::ostream &ost, const Cubes::Polygon &pg)
{
    ost << "Polygon: ";
    for (Dim p : pg)
        ost << " " << p;
    return ost << "\n";
}

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

#ifndef MEMOIZE_HASH
#define MEMOIZE_HASH 1
#endif

class Dim;
class Cells;
class PolyCube;
inline std::ostream& operator << (std::ostream&, Dim);
inline std::ostream& operator << (std::ostream&, const Cells&);
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

    struct SimpleHash
    {
        int_t operator () (const Dim &d) const
        {
            return (int_t) d.v;
        }
    };
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

    // Don't use (int_t) since PolyCube wants a symmetric hash, but we want
    // Dim{a,b} ^ Dim{c,d} != Dim{a,d} ^ Dim{c,b}. // The below hash works
    // well, i.e. we almost never get PolyCube == returns false.
    hash_t hash () const
    {
        hash_t h = 1;
        for (int i = 0; i < size (); ++i)
        {
            hash_t w0 = v[i] + 20, w = 1;
            for (int j = 0; j < i + 2; ++j)
                w *= w0;
            h = 97 + h * w;
        }
        return h;
    }

    DimIterator begin () const;
    DimIterator end () const;
};

// Iterates over x with ||x||_2 = 1.
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

struct Cells
{
    std::list<Dim> cells;

    int size () const
    {
        return cells.size ();
    }
    void add (Dim d)
    {
        const auto &end = cells.end ();
        for (auto &&p = cells.begin(); ; ++p)
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
    void erase (Dim d) // unused
    {
        const auto &end = cells.end ();
        for (auto &&p = cells.begin(); p != end; ++p)
        {
            const bool done = *p >= d;
            if (*p == d)
                cells.erase (p);
            if (done)
                break;
        }
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
    Cells operator - (Dim dim) const
    {
        Cells l (*this);
        for (Dim &d : l.cells)
            d -= dim;
        return l;
    }
};

#define BACK "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b"


struct PolyCube
{
    struct Hash
    {
        unsigned operator () (const PolyCube &pc) const
        {
            return pc.hash ();
        }
    };

    using Set = std::unordered_set<PolyCube, Hash>;
    using Corona = std::unordered_set<Dim, Dim::SimpleHash>;
    friend std::ostream& operator << (std::ostream&, const PolyCube::Corona&);

    static void merge (Set &sout, Set &sin)
    {
        sout.merge (sin);
    }
    static void merge (Set &s, PolyCube &&pc)
    {
        s.emplace (std::move (pc));
    }

    struct MuxSet
    {
        // std::mutex is not movable, so we have to use the
        // swap trick to get a vector of given size, see
        // https://stackoverflow.com/a/24170141/1556746
        std::mutex mux;
        Set set;
    };
    using Vector = std::vector<MuxSet>; // Indexed by corona size.
    const PolyCube *m_dad = nullptr;
    Dim m_cube;
#if MEMOIZE_HASH
    hash_t m_hash;
#endif

    PolyCube (const PolyCube *pc, Dim d)
        : m_dad(pc), m_cube(d)
#if MEMOIZE_HASH
        , m_hash(calc_hash ())
#endif
    {}

    Cells cubes () const
    {
        Cells c;
        for (auto p = this; p; p = p->m_dad)
            c.add (p->m_cube);
        return c;
    }
    Cells cubes_normalized () const
    {
        return cubes () - min_cube ();
    }
    Corona corona () const
    {
        const Cells &cs = cubes ();
        Corona cora;
        for (Dim d : cs.cells)
            for (Dim delta : d)
                if (! cs.contains (d + delta))
                    cora.insert (d + delta);
        return cora;
    }
    int corona_size () const
    {
        return (int) corona ().size ();
    }
    int cube_count () const
    {
        int cnt = 1;
        for (auto p = m_dad; p; p = p->m_dad)
            cnt += 1;
        return cnt;
    }
    // Symmetric in cubes and shift-invariant.
    bool operator == (const PolyCube &c) const
    {
        const auto c1 = cubes().cells;
        const auto c2 = c.cubes().cells;
        if (c1.size () != c2.size ())
            return false;
        const Dim offset = min_cube () - c.min_cube ();
        auto &&it2 = c2.begin ();
        for (Dim d1 : c1)
            if (d1 - *(it2++) != offset)
                return false;
        return true;
    }
    // Symmetric in cubes and shift-invariant.
    hash_t calc_hash () const
    {
        const Dim d = min_cube ();
        hash_t h = 0;
        for (auto p = this; p; p = p->m_dad)
            h += (p->m_cube - d).hash ();
        return h;
    }
    hash_t hash () const
    {
#if MEMOIZE_HASH
        return m_hash;
#else
        return calc_hash ();
#endif
    }
    Dim calc_min_cube () const
    {
        Dim d (m_cube);
        for (auto p = m_dad; p; p = p->m_dad)
            d = std::min (d, p->m_cube);
        return d;
    }
    Dim min_cube () const
    {
        return calc_min_cube ();
    }

    // Way 0
    void add_sprouts (Set &set) const
    {
        for (Dim d : corona())
            set.emplace (PolyCube (this, d));
    }

    // Way 4
    int add_sprouts_way4 (Vector &vms) const
    {
        int new_count = 0;
        for (Dim d : corona())
        {
            auto &&pc = PolyCube (this, d);
            MuxSet &slot = vms[pc.hash () % vms.size ()];
            slot.mux.lock ();
            const auto n = slot.set.size ();
            slot.set.emplace (pc);
            new_count += n != slot.set.size ();
            slot.mux.unlock ();
        }
        return new_count;
    }

    // Way 4
    static void add_sprouts_way4 (int n_cells, int n_slots,
                                  Vector &vset2, const Vector &vset)
    {
        Vector v (n_slots);
        vset2.swap (v); // Since resize() doesn't like std::mutex.

        // Only for printing stat.
        const int64_t n_cubes = ::cube_count (DIM, n_cells);
        std::atomic<int64_t> pc_count = 0;
        int64_t pc_show = pc_count;

#pragma omp parallel for schedule(dynamic)
        for (size_t j = 0; j < vset.size (); ++j)
        {
            for (const auto &pc : vset[j].set)
                pc_count += pc.add_sprouts_way4 (vset2);

            // Print stat.
            if (omp_get_thread_num () == 0)
            {
                const int64_t pcs = pc_count;
                if (pc_show == 0 || pcs - pc_show > 50000)
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


inline std::ostream& operator << (std::ostream &ost, const PolyCube::Corona &c)
{
    ost << "{#" << c.size ();
    for (auto d : c)
        ost << " " << d;
    return ost << " }";
}

inline std::ostream& operator << (std::ostream &ost, const PolyCube &pc)
{
    ost << "cubes: " << pc.cubes() << "\n";
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

inline std::ostream& operator << (std::ostream &ost, const Cells &c)
{
    ost << "{#" << c.size ();
    for (auto c : c.cells)
        ost << " " << c;
    return ost << " }";
}

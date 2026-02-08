// -*- c++ -*-
#include <utility> // std::move
#include <iostream>
#include <mutex>
#include <atomic>
#include <string>
#include <sstream>
#include <functional>
// Containers
#include <map>
#include <vector>
#include <unordered_set>
// C
#include <cstdint>
#include <cinttypes>
#include <cassert>
#include <climits>
#include <cstdio>
// Other
#include <omp.h>
// Own
#include "polycube-count.h"
#include "progress.h"
#include "bool-counter.h"

#include "hash.h"
#include "dim.h"
#include "cubes.h"
#include "corona.h"

BoolCounter eqne;

inline int max_possible_corona;

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

    PolyCube (const PolyCube *dad, Dim d)
        : m_cubes(dad ? & dad->m_cubes : nullptr, d)
        , m_hash (m_cubes.hash ())
    {}

    struct Iterator
    {
        CubesIterator it;
        void operator ++ () { ++it; };
        bool operator == (const Iterator &r) const { return it == r.it; }
        bool operator != (const Iterator &r) const { return it != r.it; }
        Dim operator * () const { return *it; }
    };
    Iterator begin () const  { return Iterator { m_cubes.cbegin () }; }
    Iterator end () const    { return Iterator { m_cubes.cend () }; }
    Iterator cbegin () const { return Iterator { m_cubes.cbegin () }; }
    Iterator cend () const   { return Iterator { m_cubes.cend () }; }
    using const_iterator = Iterator;

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

    struct MuxSet
    {
        // std::mutex is not movable, so we have to use the
        // swap trick to get a vector of given size, see
        // https://stackoverflow.com/a/24170141/1556746
        std::mutex mux;
        Set set;
    };
    using Vector = std::vector<MuxSet>;
private:
    Cubes m_cubes;
    hash_t m_hash = 0;
public:

    const Cubes& cubes () const
    {
        return m_cubes;
    }
    bool contains (Dim d) const
    {
        for (Dim x : *this)
            if (x == d)
                return true;
        return false;
    }
    Corona corona () const
    {
        Corona cora;
        for (Dim d : *this)
            for (Dim delta : d)
                if (! contains (d + delta))
                    cora.add (d + delta);
        return cora;
    }
    Box bounding_box () const
    {
        return m_cubes.bounding_box ();
    }
    bool has_large_corona (int max_corona) const
    {
        const Corona &cora = corona ();
        if (cora.size () > max_corona)
            return true;
        // O.b.d.A. the bounding box is oriented in a canonical way.
        const Box bb = bounding_box ();
        const Dim diam = bb.hi - bb.lo;
        for (int j = 1; j < diam.size (); ++j)
            if (diam[j] > diam[j - 1])
            return true;

        // Try some simple convexity tests.  Use an unordered_set for faster
        // accesses.  Additional size doesn't matter here.
        Corona cubs;
        for (Dim d : *this)
            cubs.add (d);
        for (Dim d : cora)
        {
            if (cubs.contains (d + Dim{1,0}) && cubs.contains (d - Dim{1,0}))
                return true;
            if (cubs.contains (d + Dim{0,1}) && cubs.contains (d - Dim{0,1}))
                return true;
        }
        return false; // Only weakly false, i.e. not necessarily convex.
    }

    bool operator == (const PolyCube &c) const
    {
        return eqne += m_cubes == c.m_cubes;
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
        for (Dim d : corona())
        {
            PolyCube pc (this, d);
            set.emplace (std::move (pc));
        }
    }

    using Filter = std::function<bool (const PolyCube&)>;

    // Way 4, 5
    int add_sprouts_way4_5 (Vector &vms, Filter filter) const
    {
        int new_count = 0;
        for (Dim d : corona())
        {
            PolyCube pc (this, d);
            if (filter && ! filter (pc))
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
                                  int progress_at, Filter filter)
    {
        Vector v (n_slots);
        vset2.swap (v); // Since resize() doesn't like std::mutex

        // Only for printing stat.
        const int64_t n_cubes = cube_count (DIM, n_cells);
        std::atomic<int64_t> pc_count = 0;
        Progress<int64_t> pro (PROGRESS_WITH_TOTAL, n_cubes, progress_at,
                               "%" PRIi64 " Cubs = %.1f%%");

#pragma omp parallel for schedule(dynamic)
        for (size_t j = 0; j < vset.size (); ++j)
        {
            for (const auto &pc : vset[j].set)
                pc_count += pc.add_sprouts_way4_5 (vset2, filter);

            // Print stat.
            if (omp_get_thread_num () == 0)
                pro.update (pc_count);
        } // parallel for
    }

    struct Poly;
    static Poly get_sprouts_poly_way5 (int N_cells, int n_cells, int n_slots,
                                       int n_parts,
                                       Vector &vset2, const Vector &vset)
    {
        if (n_cells < N_cells)
        {
            add_sprouts_way4 (n_cells, n_slots, vset2, vset, 100000, nullptr);
            return Poly (vset2);
        }
        assert (n_cells == N_cells);
        assert (n_slots >= 1);
        assert (n_parts >= 1);

        // Piecing together the poly by doing one slot at a time.
        Poly poly;
        n_slots = 2 + n_slots / n_parts;
        for (int part = 0; part < n_parts; ++part)
        {
            const Filter filter = [part, n_parts] (const PolyCube &pc)
            {
                return part == (int) (pc.hash () % n_parts);
            };
            std::cout << "part " << (1 + part) << "/" << n_parts << "\n";
            add_sprouts_way4 (n_cells, n_slots, vset2, vset, 10000, filter);
            poly += Poly (vset2);
            // This is why we are here: Purging the output each time is
            // slow but saves the memory.
            vset2.clear ();
        }
        return poly;
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

        int64_t operator () (int64_t x) const
        {
            assert (x == 1 && "todo");
            int64_t y = 0;
            for (const auto &mono : a_)
                y += mono.second;
            return y;
        }

#pragma omp declare                                                     \
    reduction (+ : Poly : omp_out += omp_in)                            \
    initializer (omp_priv = omp_orig)

#pragma omp declare                                                     \
    reduction (vecadd : std::vector<int> : vecadd (omp_out, omp_in))    \
    initializer (omp_priv = omp_orig)

        static void vecadd (std::vector<int> &a, const std::vector<int> &b)
        {
            assert (a.size () == b.size ());
            for (size_t j = 0; j < a.size (); ++j)
                a[j] += b[j];
        }

        static void init (Poly &poly, const Vector &vms)
        {
            // Way 4.
#if 0
#pragma omp parallel for schedule(dynamic,20) reduction(+: poly)
            for (size_t j = 0; j < vms.size (); ++j)
                poly += Poly (vms[j].set);
#elif 1
            std::vector<std::atomic<int>> v (1 + max_possible_corona);
            for (size_t j = 0; j < v.size (); ++j)
                v[j] = 0;
#pragma omp parallel for schedule(dynamic)
            for (size_t j = 0; j < vms.size (); ++j)
            {
                std::vector<int> w (1 + max_possible_corona, 0);
                for (const auto &pc : vms[j].set)
                {
                    const int coro = pc.corona().size ();
                    assert (coro < (int) v.size ());
                    w[coro] += 1;
                }
                for (size_t j = 0; j < v.size (); ++j)
                    if (w[j])
                        v[j] += w[j];
            }

            for (size_t j = 0; j < v.size (); ++j)
                if (v[j] != 0)
                    poly.a_[j] = v[j];
#else
            std::vector<int> v (1 + max_possible_corona, 0);
#pragma omp parallel for schedule(dynamic) reduction(vecadd: v)
            for (size_t j = 0; j < vms.size (); ++j)
            {
                std::vector<int> w (1 + max_possible_corona, 0);
                for (const auto &pc : vms[j].set)
                {
                    const int coro = pc.corona().size ();
                    assert (coro < (int) v.size ());
                    w[coro] += 1;
                }
                vecadd (v, w);
            }

            for (size_t j = 0; j < v.size (); ++j)
                if (v[j] != 0)
                    poly.a_[j] = v[j];
#endif
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
    friend std::ostream& operator << (std::ostream&, const PolyCube&);
}; // PolyCube


inline std::ostream& operator << (std::ostream &ost, const PolyCube &pc)
{
    ost << "cubes: " << pc.cubes() << "\n";
    ost << "coron: " << pc.corona() << "\n";
    return ost;
}

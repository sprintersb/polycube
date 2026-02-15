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

using Pro64 = Progress<int64_t>;
Pro64::Printer pprint64 =
    [] (std::ostringstream &ost, Pro64 &p, int64_t i) -> void
    {
        ost << i << " Cubs";
        if (p.total > 0)
            ost << " = " << 100.0 * i / p.total << "%";
    };

struct PolyCube : Cubes
{
    struct Hash
    {
        Hasher::type operator () (const PolyCube &pc) const
        {
            return pc.hash ();
        }
    };
    using Set = std::unordered_set<PolyCube, Hash>;

    static inline bool canonical_only;

    PolyCube (const PolyCube *dad, Dim d) : Cubes (dad, d) {}

#pragma omp declare reduction(merge : Set : merge (omp_out, omp_in))    \
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

public:

    void canonicalize ()
    {
        Cubes::canonicalize ();
    }
    Corona corona () const
    {
        Corona cora;
        for (Dim d : *this)
            for (Dim delta : d)
                if (! Cubes::contains (d + delta))
                    cora.add (d + delta);
        return cora;
    }
    bool has_large_corona (int max_corona) const
    {
        const Corona &cora = corona ();
        if (cora.size () > max_corona)
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

    using Filter = std::function<bool (const PolyCube&)>;

    int add_sprouts_leap_way5 (Vector &vms, int leap, Filter filter) const
    {
        assert (leap >= 2);
        int new_count = 0;
        std::vector<Set> pcs (leap);
        pcs[0].insert (*this);
        for (int i = 1; i <= leap; ++i)
        {
            for (const auto &dad : pcs[i - 1])
                for (Dim d : dad.corona())
                    if (PolyCube pc (&dad, d); i < leap)
                        pcs[i].emplace (std::move (pc));
                    else
                    {
                        if (PolyCube::canonical_only)
                            pc.canonicalize ();
                        if (! filter || filter (pc))
                        {
                            MuxSet &slot = vms[pc.hash () % vms.size ()];

                            slot.mux.lock ();
                            const auto n = slot.set.size ();
                            slot.set.emplace (std::move (pc));
                            new_count += n != slot.set.size ();
                            slot.mux.unlock ();
                        }
                    }
        }
        return new_count;
    }

    // Way 0
    void add_sprouts (Set &set) const
    {
        for (Dim d : corona())
        {
            PolyCube pc (this, d);
            if (PolyCube::canonical_only)
                pc.canonicalize ();
            set.emplace (std::move (pc));
        }
    }

    static int64_t expected_count (int n_cells)
    {
        return PolyCube::canonical_only
            ? cube_count_canon (DIM, n_cells)
            : cube_count (DIM, n_cells);
    }
    // Way 4, 5
    int add_sprouts_way4_5 (Vector &vms, Filter filter) const
    {
        int new_count = 0;
        for (Dim d : corona())
        {
            PolyCube pc (this, d);
            if (PolyCube::canonical_only)
                pc.canonicalize ();
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
    static int64_t add_sprouts_way4 (int n_cells, int n_slots, int leap,
                                     Vector &vset2, const Vector &vset,
                                     int progress_at, Filter filter,
                                     Pro64::Printer pp)
    {
        progress_at *= omp_get_max_threads ();
        Vector v (n_slots);
        vset2.swap (v); // Since resize() doesn't like std::mutex

        // Only for printing stat.
        const int64_t n_cubes = PolyCube::expected_count (n_cells);
        std::atomic<int64_t> pc_count = 0;
        if (! filter && n_cubes > 0 && leap == 0)
            std::cout << n_cubes
                      << (PolyCube::canonical_only ? " free" : " fixed")
                      << " cubs expected\n";
        Pro64 pro (n_cubes, progress_at, pp && n_cubes > 0 ? pp : pprint64);

#pragma omp parallel for schedule(dynamic)
        for (size_t j = 0; j < vset.size (); ++j)
        {
            for (const auto &pc : vset[j].set)
                pc_count += leap
                    ? pc.add_sprouts_leap_way5 (vset2, leap, filter)
                    : pc.add_sprouts_way4_5 (vset2, filter);

            // Print stat.
            if (omp_get_thread_num () == 0)
                pro.update (pc_count);
        } // parallel for
        pro.done ();
        if (! filter)
        {
            if (leap == 0)
                std::cout << pc_count
                          << (PolyCube::canonical_only ? " free" : " fixed")
                          << " cubs found\n";
            if (n_cubes > 0 && leap == 0 && pc_count != n_cubes)
                std::cout << "error: expected " << n_cubes << " != "
                          << pc_count << "\n";
        }
        return pc_count;
    }

    struct Poly;
    static Poly get_sprouts_poly_way5 (int n_cells, int n_slots,
                                       int n_parts, int leap,
                                       Vector &vset2, const Vector &vset)
    {
        if (! leap)
        {
            add_sprouts_way4 (n_cells, n_slots, 0, vset2, vset,
                              200, nullptr, nullptr);
            return Poly (vset2);
        }
        const int64_t n_cubes = PolyCube::expected_count (n_cells);
        if (n_cubes > 0)
            std::cout << n_cubes
                      << (PolyCube::canonical_only ? " free" : " fixed")
                      << " cubs expected\n";
        std::cout << "leaping " << leap << ": " << (n_cells - leap)
                  << ".." << n_cells << " with " << n_parts << " parts\n";
        std::cout.flush ();
        assert (leap >= 1);
        assert (n_slots >= 1);
        assert (n_parts >= 1);
        // leap = 1 isn't real leap-frogging.
        leap = leap == 1 ? 0 : leap;

        // Piecing together the poly by doing one slot at a time.
        Poly poly;
        n_slots = 2 + n_slots / n_parts;
        Pro64::Printer pbar =
            [n_parts] (std::ostringstream &ost, Pro64 &p, int64_t i) -> void
            {
                p.print_bar (ost, 79, (double) n_parts * i / p.total);
            };
        int64_t pc_count = 0;
        for (int part = 0; part < n_parts; ++part)
        {
            const Filter filter = [part, n_parts] (const PolyCube &pc)
            {
                return part == (int) (pc.hash () % n_parts);
            };
            std::cout << "part " << (1 + part) << "/" << n_parts << "\n";
            pc_count += add_sprouts_way4 (n_cells, n_slots, leap, vset2, vset,
                                          1 + 200 / n_parts, filter, pbar);
            poly += Poly (vset2);
            // This is why we are here: Purging the output each time is
            // slow but saves the memory.
            vset2.clear ();
        }
        std::cout << pc_count
                  << (PolyCube::canonical_only ? " free" : " fixed")
                  << " cubs found\n";
        if (n_cubes > 0 && pc_count != n_cubes)
            std::cout << "error: expected " << n_cubes << " != "
                      << pc_count << "\n";
        return poly;
    }

    // Univariate polynomial over Z in sparse representation.
    struct Poly
    {
        using T = int64_t;
        using value_type = T;
        using coefficient_type = T;
        using exponent_type = int;
        using monomial_type = std::map<int,T>::value_type;
        enum { POLY_TEX, POLY_LIST };
        std::map<int,T> a_;

        Poly () {}

        Poly (const monomial_type &mono)
        {
            a_[mono.first] = mono.second;
        }
        // Way 0: 100% sequential.
        Poly (const Set &set)
        {
            for (const auto &pc : set)
            {
                int multi = PolyCube::canonical_only ? pc.multiplicity() : 1;
                const int coro = pc.corona().size ();
                const auto &monome = a_.find (coro);
                if (monome == a_.end ())
                    a_[coro] = multi;
                else
                    monome->second += multi;
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
                const T coef = q_mono.second;
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
    reduction (vecadd : std::vector<T> : vecadd (omp_out, omp_in))      \
    initializer (omp_priv = omp_orig)

        static void vecadd (std::vector<T> &a, const std::vector<T> &b)
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
            std::vector<std::atomic<T>> v (1 + max_possible_corona);
            for (size_t j = 0; j < v.size (); ++j)
                v[j] = 0;
#pragma omp parallel for schedule(dynamic)
            for (size_t j = 0; j < vms.size (); ++j)
            {
                std::vector<T> w (1 + max_possible_corona, 0);
                for (const auto &pc : vms[j].set)
                {
                    const int coro = pc.corona().size ();
                    assert (coro < (int) v.size ());
                    w[coro] += PolyCube::canonical_only ? pc.multiplicity() : 1;
                }
                for (size_t j = 0; j < v.size (); ++j)
                    if (w[j])
                        v[j] += w[j];
            }

            for (size_t j = 0; j < v.size (); ++j)
                if (v[j] != 0)
                    poly.a_[j] = v[j];
#else
            std::vector<T> v (1 + max_possible_corona, 0);
#pragma omp parallel for schedule(dynamic) reduction(vecadd: v)
            for (size_t j = 0; j < vms.size (); ++j)
            {
                std::vector<T> w (1 + max_possible_corona, 0);
                for (const auto &pc : vms[j].set)
                {
                    const int coro = pc.corona().size ();
                    assert (coro < (int) v.size ());
                    w[coro] += PolyCube::canonical_only ? pc.multiplicity() : 1;
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
                        printf ("%" PRIi64, m.second);
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
                    printf ("%" PRIi64 ",%d", m.second, m.first);
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
    ost << "cubes: " << static_cast<Cubes> (pc) << "\n";
    ost << "coron: " << pc.corona() << "\n";
    return ost;
}

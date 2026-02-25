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
#include <climits>
#include <cstdio>
// Other
#include <omp.h>
// Own
#include "oeis.h"
#include "progress.h"
#include "diagnostic.h"
#include "debug.h"
#include "hash.h"
#include "cubes.h"
#include "corona.h"

inline int max_possible_corona;

using Pro64 = Progress<int64_t>;

inline Pro64::Printer pprint64 =
    [] (std::ostringstream &ost, Pro64 &p, int64_t i) -> void
    {
        ost << i << " Cubs";
        if (p.total > 0)
            ost << " = " << 100.0 * i / p.total << "%";
    };

inline Pro64::Printer pprint64_int =
    [] (std::ostringstream &ost, Pro64 &p, int64_t i) -> void
    {
        p.print_bar (ost, 79, (double) i / p.total);
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

    PolyCube (const std::initializer_list<Dim> &il) : Cubes (il) {}
    PolyCube (const PolyCube &dad, Dim d) : Cubes (dad, d) {}
    PolyCube (const Cubes &c) : Cubes (c) {}

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
    using Filter = std::function<bool (const PolyCube&)>;

#include "polynomial.h"

    struct Flavour
    {
        Flavour () = delete;
        std::string name;
        Flavour (const char *name) : name(std::string (name)) {}
        std::optional<int64_t> count;
        Set set;
        // Algorithm works on many hash slots to avoid thread clogging.
        Vector vec;
        bool has_set () const
        {
            return ! set.empty ();
        }
        bool has_vec () const
        {
            for (const auto &ms : vec)
                if (! ms.set.empty ())
                    return true;
            return false;
        }
        bool is_free () const { return name == "free"; }
        bool is_fixed () const { return name == "fixed"; }
        friend std::ostream& operator << (std::ostream&, const Flavour&);
    };
    struct Have
    {
        Flavour fixed { "fixed" }, free { "free" };
        Flavour symm { "symm" }, asymm { "asymm" };
        //Have
        Flavour& operator () ()
        {
            return canonical_only ? free : fixed;
        }
        const Flavour& operator () () const
        {
            return canonical_only ? free : fixed;
        }
        std::optional<Poly> corona_polynomial;
        std::optional<Cubes> smallest_corona;
        int64_t get_fixed_count () const;
        int64_t get_free_count () const;
    };

    struct Want
    {
        struct
        {
            bool count, cubes;
        } free, fixed;
        bool corona_polynomial;
        bool smallest_corona;
        int way;
        int leap;
        int n_parts;
        int n_cells;
        int n_cells_final;
        bool progress;
    };
    static inline Want want;

public:

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

    int sprout_leap_way5 (Vector &vms, int leap, const Filter &filter) const
    {
        Assert (leap >= 2, "leap %d is no real leapfrogging", leap);
        int new_count = 0;
        std::vector<Set> pcs (leap);
        pcs[0].insert (*this);
        for (int i = 1; i <= leap; ++i)
        {
            for (const auto &dad : pcs[i - 1])
                for (Dim d : dad.corona())
                    if (PolyCube pc (dad, d); i < leap)
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
            PolyCube pc (*this, d);
            if (PolyCube::canonical_only)
                pc.canonicalize ();
            set.emplace (std::move (pc));
        }
    }

    static int64_t expected_count (int n_cells)
    {
        return PolyCube::canonical_only
            ? oeis::cubes_free (DIM, n_cells)
            : oeis::cubes_fixed (DIM, n_cells);
    }
    // Way 4, 5
    int sprout_way4_5 (Vector &vms, const Filter &filter) const
    {
        int new_count = 0;
        for (Dim d : corona())
        {
            PolyCube pc (*this, d);
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
    static int64_t sprout_way4 (int n_slots, int progress_at,
                                Have &hout, const Have &hin,
                                const Filter &filter, Pro64::Printer pp)
    {
        progress_at *= omp_get_max_threads ();
        if ((int) hout().vec.size () < n_slots)
        {
            Vector v (n_slots);
            hout().vec.swap (v); // Since resize() doesn't like std::mutex
        }

        // Only for printing stat.
        const int64_t n_cubes = PolyCube::expected_count (want.n_cells);
        std::atomic<int64_t> pc_count = 0;
        if (! filter && n_cubes > 0 && want.leap == 0)
            std::cout << n_cubes << " " << hout() << " cubs expected\n";
        Pro64 pro (n_cubes, progress_at, pp && n_cubes > 0 ? pp : pprint64);

#pragma omp parallel for schedule(dynamic)
        for (size_t j = 0; j < hin().vec.size (); ++j)
        {
            for (const auto &pc : hin().vec[j].set)
                pc_count += want.leap >= 2
                    ? pc.sprout_leap_way5 (hout().vec, want.leap, filter)
                    : pc.sprout_way4_5 (hout().vec, filter);
            // Print stat.
            if (want.progress && omp_get_thread_num () == 0)
                    pro.update (pc_count);
        } // parallel for
        pro.done ();
        if (! filter)
        {
            if (want.leap == 0)
            {
                std::cout << pc_count << " " << hout() << " cubs found\n";
                hout().count = pc_count;
            }
            if (n_cubes > 0 && want.leap == 0 && pc_count != n_cubes)
                std::cout << "error: expected " << n_cubes << " != "
                          << pc_count << "\n";
        }
        return pc_count;
    }

    struct Poly;
    static void sprout_way5 (int n_slots, int progress_at,
                             Have &hout, const Have &hin)
    {
        if (! want.leap)
        {
            sprout_way4 (n_slots, progress_at, hout, hin, nullptr, nullptr);
            if (want.corona_polynomial)
                hout.corona_polynomial = Poly (hout().vec);
            return;
        }
        const int64_t n_cubes = PolyCube::expected_count (want.n_cells);
        int leap = want.leap;
        if (n_cubes > 0)
            std::cout << n_cubes << " " << hout() << " cubs expected\n";
        std::cout << "leaping " << leap << ": " << (want.n_cells - leap)
                  << ".." << want.n_cells << " with " << want.n_parts
                  << " parts\n";
        std::cout.flush ();
        Assert (want.n_parts >= 1, "cannot leap with %d parts", want.n_parts);
        // leap = 1 isn't real leap-frogging.
        leap = leap == 1 ? 0 : leap;

        if (want.free.cubes || want.fixed.cubes)
            error ("leaping can't produce cubes, only counts and polynomial");

        // Piecing together the poly by doing one slot at a time.
        if (want.corona_polynomial)
            hout.corona_polynomial = Poly ();
        auto &poly = hout.corona_polynomial;
        auto &count = hout().count;
        count = 0;

        n_slots = 2 + n_slots / want.n_parts;
        Pro64::Printer pbar =
            [] (std::ostringstream &ost, Pro64 &p, int64_t i) -> void
            {
                p.print_bar (ost, 79, (double) want.n_parts * i / p.total);
            };

        if (want.fixed.count && hout().is_free ())
        {
            if (Cubes::can_canonicalize)
                hout.fixed.count = 0;
            else
                error ("can't canonicalize in dim %d", DIM);
        }

        for (int part = 0; part < want.n_parts; ++part)
        {
            const Filter filter = [part] (const PolyCube &pc)
            {
                return part == (int) (pc.hash () % want.n_parts);
            };
            std::cout << "part " << (1 + part) << "/" << want.n_parts << "\n";
            *count += sprout_way4 (n_slots, 1 + progress_at / want.n_parts,
                                   hout, hin, filter, pbar);
            if (poly)
                *poly += Poly (hout().vec);
            if (want.fixed.count && hout().is_free ())
            {
                std::cout << "free -> fixed\n";
                int64_t cnt = 0;
                std::atomic<int> i_done = 0;
                Pro64 p (hout().vec.size(), omp_get_max_threads() / 2,
                         pprint64_int);
#pragma omp parallel for schedule(dynamic) reduction(+ : cnt)
                for (size_t i = 0; i < hout().vec.size (); ++i)
                {
                    uint64_t c = 0;
                    for (const auto &pc : hout().vec[i].set)
                        c += pc.multiplicity ();
                    cnt += c;
                    i_done += 1;
                    if (want.progress && omp_get_thread_num () == 0)
                        p.update (i_done);
                }
                * hout.fixed.count += cnt;
            }
            // This is why we are here: Purging the output each time is
            // slow but saves the memory.
#pragma omp parallel for schedule(dynamic)
            for (size_t i = 0; i < hout().vec.size (); ++i)
                hout().vec[i].set.clear ();
        }
        hout().vec.clear ();

        std::cout << *count << " " << hout() << " cubs found\n";
        if (n_cubes > 0 && *count != n_cubes)
            std::cout << "error: expected " << n_cubes << " != "
                      << *count << "\n";
    }
    friend std::ostream& operator << (std::ostream&, const PolyCube&);
}; // PolyCube


inline std::ostream& operator << (std::ostream &ost, const PolyCube &pc)
{
    ost << "cubes: " << static_cast<Cubes> (pc) << "\n";
    ost << "coron: " << pc.corona() << "\n";
    return ost;
}

inline std::ostream& operator << (std::ostream &ost,
                                  const PolyCube::Flavour &flavour)
{
    return ost << flavour.name;
}

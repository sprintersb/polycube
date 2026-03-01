// -*- c++ -*-
#include <utility> // std::move
#include <iostream>
#include <mutex>
#include <atomic>
#include <string>
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
#include "progress.h"
#include "hash.h"
#include "cubes.h"
#include "corona.h"

using Pro64 = Progress<int64_t>;

inline int max_possible_corona;

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

    struct MuxSet
    {
        // std::mutex is not movable, so we have to use the
        // swap trick to get a vector of given size, see
        // https://stackoverflow.com/a/24170141/1556746
        std::mutex mux;
        Set set;
        bool insert (const PolyCube&);
    };
    using VectorBase = std::vector<MuxSet>;
    struct Vector : VectorBase
    {
        Vector () : VectorBase () {}
        Vector (size_t s) : VectorBase (s) {}
        // Since resize() doesn't like std::mutex.
        void resize (size_t sz)
        {
            if (sz != size ())
            {
                Vector v (sz);
                VectorBase::swap (v);
            }
        }
    };
    struct Have;
    struct Flavour
    {
        Flavour () = delete;
        const std::string name;
        const bool is_fixed, is_symm0, is_symm1, is_symm2;
        Cubes::Symmetry symm;
        Flavour (const char *name) : name(std::string (name))
            , is_fixed(! strcmp (name, "fixed"))
            , is_symm0(! strcmp (name, "symm0"))
            , is_symm1(! strcmp (name, "symm1"))
            , is_symm2(! strcmp (name, "symm2"))
        {
            const int id = name[4] - '0';
            symm = id >= 0 && id <= 2 ? id : 0;
        }
        std::optional<int64_t> count;
        // Algorithm works on many hash slots to avoid thread clogging.
        Vector vec;
        bool has_vec () const
        {
            for (const auto &ms : vec)
                if (! ms.set.empty ())
                    return true;
            return false;
        }
        friend std::ostream& operator << (std::ostream&, const Flavour&);
    };
    using VOut = std::vector<Flavour*>;
    using VIn = std::vector<const Flavour*>;
    using Filter = std::function<bool (const PolyCube&)>;

    // Univariate polynomial over Z in sparse representation.
    struct Poly
    {
        using T = int64_t;
        using value_type = T;
        using coefficient_type = T;
        using exponent_type = int;
        using monomial_type = std::map<int,T>::value_type;
        enum { POLY_TEX, POLY_LIST };
        std::map<int, T> a_;

        Poly () {}
        Poly (const monomial_type &mono)
        {
            a_[mono.first] = mono.second;
        }
        Poly (const Have &have)
        {
            for (auto f : have.v ())
                *this += Poly (f);
        }
        Poly (const Flavour *fl)
        {
            Poly::init (*this, fl);
        }
        static void init (Poly&, const Flavour*);
        Poly& operator += (const Poly&);
        int64_t operator () (int64_t) const;
        void print (int n, int style, const char *var = "q") const;
    }; // Poly

    struct Have
    {
        Flavour fixed { "fixed" };
        Flavour symm0 { "symm0" }, symm1 { "symm1" }, symm2 { "symmm2" };
        std::optional<Poly> corona_polynomial;
        std::optional<Cubes> smallest_corona;
        std::optional<int64_t> free_count;

        void init (const PolyCube&, int n_slots);
        int64_t get_fixed_count () const;
        int64_t get_free_count () const;
        VOut v ()
        {
            return canonical_only
                ? VOut { &symm0, &symm1, &symm2 }
                : VOut { &fixed };
        }
        VIn v () const
        {
            return canonical_only
                ? VIn { &symm0, &symm1, &symm2 }
                : VIn { &fixed };
        }
        int64_t sprout_way4 (int n_slots, int progress_at, Have &hout,
                             const Filter &filter, Pro64::Printer pp) const;
        void sprout_way5 (int n_slots, int progress_at, Have &hout) const;
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

    static Set find_min_corona (const Have &have)
    {
        Set set;
        std::atomic<int> corona_size = INT_MAX;
        for (auto f : have.v ())
#pragma omp parallel for schedule(dynamic) // reduction(min: corona_size)
            for (size_t j = 0; j < f->vec.size (); ++j)
                for (const auto &pc : f->vec[j].set)
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

    static int64_t expected_count (int n_cells);

    bool insert (VOut &vout, const Filter &filter);
    int sprout_leap (VOut &vout, int leap, const Filter &filter) const;
    int sprout (VOut &vout, const Filter &filter) const;

    friend std::ostream& operator << (std::ostream&, const PolyCube&);
}; // PolyCube

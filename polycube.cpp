#include <set>
#include <omp.h>
// Own
#include "polycube.h"
#include "oeis.h"
#include "diagnostic.h"
#include "debug.h"

/////////////////////////////////////////////////////////////////////////////
// Sprouting

int64_t PolyCube::expected_count (int n_cells)
{
    return PolyCube::canonical_only
        ? oeis::cubes_free (DIM, n_cells)
        : oeis::cubes_fixed (DIM, n_cells);
}

void PolyCube::Want::announce_expectations () const
{
    const int64_t n_cubes = PolyCube::expected_count (want.n_cells);
    const char *kind = PolyCube::canonical_only ? "free" : "fixed";
    if (n_cubes > 0)
        std::cout << n_cubes << " " << kind << " cubs expected\n";
}

void PolyCube::Have::show_result () const
{
    const int64_t n_cubes = PolyCube::expected_count (want.n_cells);
    const char *kind = PolyCube::canonical_only ? "free" : "fixed";
    int64_t ccount = PolyCube::canonical_only
        ? *count.free : *fixed.count;

    std::cout << ccount << " " << kind << " cubs found\n";
    if (n_cubes > 0 && ccount != n_cubes)
        std::cout << "error: expected " << n_cubes << " != "
                  << ccount << "\n";

    if (canonical_only)
    {
        auto s1 = oeis::cubes ("asymm", DIM);
        auto a1 = s1[want.n_cells];
        auto c1 = *count.asymm;
        std::cout << s1.name << " = " << a1 << " = " << c1 << "\n";
        auto s2 = oeis::cubes ("symm", DIM);
        auto a2 = s2[want.n_cells];
        auto c2 = *count.symm;
        std::cout << s2.name << " = " << a2 << " = " << c2 << "\n";
    }
}

inline Pro64::Printer pprint64 =
    [] (std::ostringstream &ost, Pro64 &p, int64_t i) -> void
    {
        if (p.total > 0)
            p.print_bar (ost, 79, (double) i / p.total);
    };

inline Pro64::Printer pprint64_int =
    [] (std::ostringstream &ost, Pro64 &p, int64_t i) -> void
    {
        p.print_bar (ost, 79, (double) i / p.total);
    };

inline bool PolyCube::MuxSet::insert (const PolyCube &pc)
{
    mux.lock ();
    const auto n = set.size ();
    set.emplace (pc);
    const bool added = n != set.size ();
    mux.unlock ();
    return added;
}

inline bool PolyCube::insert (VOut &vout, const Filter &filter)
{
    Symmetry symm = 0;
    if (PolyCube::canonical_only)
        canonicalize (symm);
    if (! filter || filter (*this))
    {
        Vector &vms = vout[symm]->vec;
        MuxSet &slot = vms[hash () % vms.size ()];
        return slot.insert (*this);
    }
    return false;
}

// Add filtered sprouts to vout and return the number of adds.
int PolyCube::sprout_leap (VOut &vout, int leap, const Filter &filter) const
{
    Assert (leap >= 2, "leap %d is no real leapfrogging", leap);
    int new_count = 0;
    std::vector<Set> pcs (leap);
    pcs[0].insert (*this);
    for (int i = 1; i <= leap; ++i)
        for (const auto &dad : pcs[i - 1])
            for (Dim d : dad.corona())
                if (PolyCube pc (dad, d); i < leap)
                    pcs[i].emplace (std::move (pc));
                else
                    new_count += pc.insert (vout, filter);
    return new_count;
}

// Add filtered sprouts to vout and return the number of adds.
int PolyCube::sprout (VOut &vout, const Filter &filter) const
{
    int new_count = 0;
    for (Dim d : corona())
        new_count += PolyCube (*this, d).insert (vout, filter);
    return new_count;
}

int64_t PolyCube::Have::sprout (int n_slots, const Have &hin,
                                const Filter &filter, Pro64::Printer pp)
{
    const auto vin = hin.v();
    auto vout = v();
    for (auto f : vout)
        f->vec.resize (n_slots);

    // Progress: "Integrate" over the input to get exact mileage.
    std::vector<int64_t> dd;
    int64_t n_in = 0;
    for (auto f : vin)
        for (auto &ms : f->vec)
            dd.push_back (n_in += ms.set.size ());
    Pro64 pro (n_in, 1, pp ? pp : pprint64);

    std::atomic<int64_t> pc_count = 0;
    int64_t j_off = 0;
    for (auto f : vin)
    {
#pragma omp parallel for schedule(dynamic)
        for (size_t j = 0; j < f->vec.size (); ++j)
        {
            for (const auto &pc : f->vec[j].set)
                pc_count += want.leap >= 2
                    ? pc.sprout_leap (vout, want.leap, filter)
                    : pc.sprout (vout, filter);
            // Print stat.
            if (want.progress && omp_get_thread_num () == 0)
                pro.update (dd[j + j_off]);
        } // parallel for
        j_off += f->vec.size ();
    }
    pro.done ();
    for (auto f : vout)
    {
        f->count = 0;
        for (auto &ms : f->vec)
            *f->count += ms.set.size ();
    }

    return pc_count;
}


void PolyCube::Have::sprout_way5 (int n_slots, const Have &hin)
{
    int64_t ccount = 0;
    auto &poly = corona_polynomial;

    if (! want.leap)
    {
        ccount = sprout (n_slots, hin, nullptr, nullptr);
        if (want.corona_polynomial)
            poly = Poly (*this);
        if (!! symm0.count && !! symm1.count && !! symm2.count)
        {
            count.symm = *symm2.count;
            count.asymm = *symm0.count + *symm1.count;
            count.free = *count.asymm + *count.symm;
        }
    }
    else
    {
        const auto vin = hin.v();
        auto vout = v();
        int leap = want.leap;
        std::cout << "want:" << (want.free.count ? " free.count" : "")
                  << (want.fixed.count ? " fixed.count" : "")
                  << (want.corona_polynomial ? " polynomial" : "") << "\n";

        std::cout << "leaping " << leap << ": " << (want.n_cells - leap)
                  << ".." << want.n_cells << " with " << want.n_parts
                  << " parts\n";
        std::cout.flush ();
        Assert (want.n_parts >= 1, "cannot leap with %d parts", want.n_parts);
        // leap = 1 isn't real leap-frogging.
        leap = leap == 1 ? 0 : leap;

        if (want.free.cubes || want.fixed.cubes)
            error ("leaping can't produce cubes, only counts and polynomial");
        if (want.free.count && ! Cubes::can_canonicalize)
            error ("can't canonicalize in dim %d", DIM);
        if (! want.free.count && ! want.fixed.count && ! want.corona_polynomial)
            error ("you better want something");

        // Piecing together the poly by doing one slot at a time.
        auto &poly = corona_polynomial;
        if (want.corona_polynomial)
            poly = Poly ();
        if (want.fixed.count)
            fixed.count = 0;
        if (want.free.count)
        {
            symm0.count = symm1.count = symm2.count = 0;
            count.symm = count.asymm = 0;
        }

        n_slots = 2 + n_slots / want.n_parts;

        for (int part = 0; part < want.n_parts; ++part)
        {
            const Pro64::Printer pbar = [part] (auto &ost, Pro64 &p, int64_t i)
            {
                ost << (1 + part) << "/" << want.n_parts << " ";
                p.print_bar (ost, 79, (double) i / p.total);
            };
            const Filter filter = [part] (const PolyCube &pc)
            {
                return part == (int) (pc.hash () % want.n_parts);
            };

            ccount += sprout (n_slots, hin, filter, pbar);

            if (PolyCube::canonical_only)
            {
                *count.asymm += *symm0.count + *symm1.count;
                *count.symm += *symm2.count;
            }
            if (poly)
                *poly += Poly (*this);
            if (want.fixed.count && PolyCube::canonical_only && !poly)
            {
                std::cout << "free -> fixed\n";
                std::atomic<int> i_done = 0;
                Pro64 p (3 * n_slots, omp_get_max_threads() / 2, pprint64_int);
                for (auto f : vout)
                {
                    int64_t cnt = 0;
#pragma omp parallel for schedule(dynamic) reduction(+ : cnt)
                    for (size_t i = 0; i < f->vec.size (); ++i)
                    {
                        uint64_t c = 0;
                        for (const auto &pc : f->vec[i].set)
                            c += pc.multiplicity (f->symm);
                        cnt += c;
                        i_done += 1;
                        if (want.progress && omp_get_thread_num () == 0)
                            p.update (i_done);
                    }
                    *fixed.count += cnt;
                }
            }
            // This is why we are here: Purging the output each time is
            // slow but saves the memory.
            for (auto f : vout)
#pragma omp parallel for schedule(dynamic)
                for (size_t i = 0; i < f->vec.size (); ++i)
                    f->vec[i].set.clear ();
        } // for n_parts.
        for (auto f : vout)
            f->vec.clear ();
    }

    if (PolyCube::canonical_only)
        count.free = ccount;
    else
        fixed.count = ccount;
    if (want.fixed.count && PolyCube::canonical_only && poly)
        fixed.count = (*poly) (1);
}


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


/////////////////////////////////////////////////////////////////////////////
// PolyCube::Have

void PolyCube::Have::init (const PolyCube &pc, int n_slots)
{
    const Cubes::Symmetry symm = PolyCube::canonical_only ? 2 : 0;
    VOut vout = v ();
    for (auto f : vout)
    {
        f->vec.resize (n_slots);
        f->count = 0;
    }
    Vector &vms = vout[symm]->vec;
    vms[pc.hash () % n_slots].set.emplace (pc);
    vout[symm]->count = 1;
}

int64_t PolyCube::Have::get_fixed_count () const
{
    int64_t cnt = 0;

    if (fixed.count)
        return *fixed.count;
    else if (fixed.has_vec ())
        for (const auto &ms : fixed.vec)
            cnt += ms.set.size ();
    else if (corona_polynomial)
        cnt = (*corona_polynomial) (1);
    else if (Cubes::can_canonicalize && symm2.has_vec ())
    {
        for (auto f : v ())
        {
            int64_t c2 = 0;
#pragma omp parallel for schedule(dynamic) reduction(+ : c2)
            for (size_t i = 0; i < f->vec.size (); ++i)
            {
                int64_t c = 0;
                for (const auto &pc : f->vec[i].set)
                    c += pc.multiplicity (f->symm);
                c2 += c;
            }
            cnt += c2;
        }
    }
    else
        warning ("don't know how to obtain fixed cubes count");

    return cnt ? cnt : -1;
}

int64_t PolyCube::Have::get_free_count () const
{
    int64_t cnt = 0;

    if (count.free)
        return *count.free;
    else if (symm2.count)
        return *symm0.count + *symm1.count + *symm2.count;
    else if (symm2.has_vec ())
    {
        for (const auto &ms : symm0.vec)
            cnt += ms.set.size ();
        for (const auto &ms : symm1.vec)
            cnt += ms.set.size ();
        for (const auto &ms : symm2.vec)
            cnt += ms.set.size ();
    }
    else if (Cubes::can_canonicalize && fixed.has_vec ())
#pragma omp parallel for schedule(dynamic) reduction(+ : cnt)
        for (size_t i = 0; i < fixed.vec.size (); ++i)
        {
            int64_t c = 0;
            for (const auto &pc : fixed.vec[i].set)
                c += pc.is_canonical ();
            cnt += c;
        }
    else
        warning ("don't know how to obtain free cubes count");

    return cnt ? cnt : -1;
}

/////////////////////////////////////////////////////////////////////////////
// Polynomial

using T = PolyCube::Poly::value_type;

static void vecadd (std::vector<T> &a, const std::vector<T> &b)
{
    Assert (a.size () == b.size (), "vector add of different sizes");
    for (size_t j = 0; j < a.size (); ++j)
        a[j] += b[j];
}

#pragma omp declare                                     \
    reduction (+ : PolyCube::Poly : omp_out += omp_in)  \
    initializer (omp_priv = omp_orig)

#pragma omp declare                                                 \
    reduction (vecadd : std::vector<T> : vecadd (omp_out, omp_in))  \
    initializer (omp_priv = omp_orig)

void PolyCube::Poly::init (Poly &poly, const Flavour *f)
    {
        std::vector<std::atomic<T>> v (1 + max_possible_corona);
        for (size_t j = 0; j < v.size (); ++j)
            v[j] = 0;
#pragma omp parallel for schedule(dynamic)
        for (size_t j = 0; j < f->vec.size (); ++j)
        {
            std::vector<T> w (1 + max_possible_corona, 0);
            for (const auto &pc : f->vec[j].set)
            {
                const int coro = pc.corona().size ();
                Assert (coro < (int) v.size (), "impossible corona size");
                w[coro] += PolyCube::canonical_only
                    ? pc.multiplicity (f->symm)
                    : 1;
            }
            for (size_t j = 0; j < v.size (); ++j)
                if (w[j])
                    v[j] += w[j];
        }

        for (size_t j = 0; j < v.size (); ++j)
            if (v[j] != 0)
                poly.a_[j] = v[j];
    }

auto PolyCube::Poly::operator += (const Poly &q) -> Poly&
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

int64_t PolyCube::Poly::operator () (int64_t x) const
{
    (void) x;
    Assert (x == 1, "todo: evaluate Poly at x != 1");
    int64_t y = 0;
    for (const auto &mono : a_)
        y += mono.second;
    return y;
}


void PolyCube::Poly::print (int n, int style, const char *var) const
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

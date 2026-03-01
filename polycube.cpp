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
    bool symm = 0;
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

int64_t PolyCube::Have::sprout_way4 (int n_slots, int progress_at, Have &hout,
                                     const Filter &filter, Pro64::Printer pp)
    const
{
    progress_at *= omp_get_max_threads ();
    const auto vin = v();
    auto vout = hout.v();
    for (auto f : vout)
        f->vec.resize (n_slots);

    // Only for printing stat.
    const int64_t n_cubes = PolyCube::expected_count (want.n_cells);
    const char *kind = PolyCube::canonical_only ? "free" : "fixed";
    if (! filter && n_cubes > 0 && ! want.leap)
        std::cout << n_cubes << " " << kind << " cubs expected\n";
    Pro64 pro (n_cubes, progress_at, pp && n_cubes > 0 ? pp : pprint64);

    std::atomic<int64_t> pc_count = 0;
    for (auto f : vin)
#pragma omp parallel for schedule(dynamic)
        for (size_t j = 0; j < f->vec.size (); ++j)
        {
            for (const auto &pc : f->vec[j].set)
                pc_count += want.leap >= 2
                    ? pc.sprout_leap (vout, want.leap, filter)
                    : pc.sprout (vout, filter);
            // Print stat.
            if (want.progress && omp_get_thread_num () == 0)
                pro.update (pc_count);
        } // parallel for
    pro.done ();
    for (auto f : vout)
    {
        f->count = 0;
        for (auto &ms : f->vec)
            *f->count += ms.set.size ();
    }

    if (! filter && ! want.leap)
    {
        std::cout << pc_count << " " << kind << " cubs found\n";
        if (n_cubes > 0 && pc_count != n_cubes)
            std::cout << "error: expected " << n_cubes << " != "
                      << pc_count << "\n";
        for (auto f : vout)
        {
            auto s = oeis::cubes (f->name, DIM);
            auto a = s[want.n_cells];
            std::cout << s.name << " = " << a << " = " << *f->count << "\n";
        }
    }
    return pc_count;
}


void PolyCube::Have::sprout_way5 (int n_slots, int progress_at,
                                  Have &hout) const
{
    if (! want.leap)
    {
        sprout_way4 (n_slots, progress_at, hout, nullptr, nullptr);
        if (want.corona_polynomial)
            hout.corona_polynomial = Poly (hout);
        return;
    }
    const auto vin = v();
    auto vout = hout.v();
    const int64_t n_cubes = PolyCube::expected_count (want.n_cells);
    int leap = want.leap;
    const char *kind = PolyCube::canonical_only ? "free" : "fixed";
    if (n_cubes > 0)
        std::cout << n_cubes << " " << kind << " cubs expected\n";
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
    auto &poly = hout.corona_polynomial;
    if (want.corona_polynomial)
        poly = Poly ();
    if (want.fixed.count)
        hout.fixed.count = 0;
    if (want.free.count)
        hout.symm.count = hout.asymm.count = 0;

    n_slots = 2 + n_slots / want.n_parts;
    Pro64::Printer pbar =
        [] (std::ostringstream &ost, Pro64 &p, int64_t i) -> void
        {
            p.print_bar (ost, 79, (double) want.n_parts * i / p.total);
        };

    int64_t count = 0;
    for (int part = 0; part < want.n_parts; ++part)
    {
        const Filter filter = [part] (const PolyCube &pc)
        {
            return part == (int) (pc.hash () % want.n_parts);
        };
        std::cout << "part " << (1 + part) << "/" << want.n_parts << "\n";
        count += sprout_way4 (n_slots, 1 + progress_at / want.n_parts,
                              hout, filter, pbar);
        if (poly)
            *poly += Poly (hout);
        if (want.fixed.count && PolyCube::canonical_only && !poly)
        {
            std::cout << "free -> fixed\n";
            std::atomic<int> i_done = 0;
            Pro64 p (2 * n_slots, omp_get_max_threads() / 2, pprint64_int);
            for (auto f : vout)
            {
                int64_t cnt = 0;
#pragma omp parallel for schedule(dynamic) reduction(+ : cnt)
                for (size_t i = 0; i < f->vec.size (); ++i)
                {
                    uint64_t c = 0;
                    for (const auto &pc : f->vec[i].set)
                        c += pc.multiplicity (f->is_symm);
                    cnt += c;
                    i_done += 1;
                    if (want.progress && omp_get_thread_num () == 0)
                        p.update (i_done);
                }
                *hout.fixed.count += cnt;
            }
        }
        // This is why we are here: Purging the output each time is
        // slow but saves the memory.
        for (auto f : vout)
#pragma omp parallel for schedule(dynamic)
            for (size_t i = 0; i < f->vec.size (); ++i)
                f->vec[i].set.clear ();
    }
    for (auto f : vout)
        f->vec.clear ();

    if (PolyCube::canonical_only)
        hout.free_count = count;
    else
        hout.fixed.count = count;
    if (want.fixed.count && PolyCube::canonical_only && poly)
        hout.fixed.count = (*poly) (1);

    std::cout << count << " " << kind << " cubs found\n";
    if (n_cubes > 0 && count != n_cubes)
        std::cout << "error: expected " << n_cubes << " != "
                  << count << "\n";
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
    const bool symm = PolyCube::canonical_only;
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
    else if (Cubes::can_canonicalize && symm.has_vec ())
    {
        for (auto f : v ())
        {
            int64_t c2 = 0;
#pragma omp parallel for schedule(dynamic) reduction(+ : c2)
            for (size_t i = 0; i < f->vec.size (); ++i)
            {
                int64_t c = 0;
                for (const auto &pc : f->vec[i].set)
                    c += pc.multiplicity (f->is_symm);
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

    if (free_count)
        return *free_count;
    else if (symm.count)
        return *symm.count + *asymm.count;
    else if (symm.has_vec ())
    {
        for (const auto &ms : symm.vec)
            cnt += ms.set.size ();
        for (const auto &ms : asymm.vec)
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
                    ? pc.multiplicity (f->is_symm)
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

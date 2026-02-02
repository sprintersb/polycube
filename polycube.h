// -*- c++ -*-
#include <array>
#include <list>
#include <map>
#include <vector>
#include <set>
#include <unordered_set>
#include <utility> // std::move
#include <iostream>
#include <iterator>
#include <mutex>

#include <cstdint>
#include <cassert>
#include <omp.h>

inline int64_t cube_count (int dim, int n_cells)
{
    static const std::array<std::vector<int64_t>, 4 + 1> cube_counts =
    {
        {
            {},
            {},
            // n = 2: https://oeis.org/A001168
            { 1, 1, 2, 6, 19, 63, 216, 760, 2725, 9910, 36446, 135268, 505861,
              1903890, 7204874, 27394666, 104592937, 400795844, 1540820542,
              5940738676, 22964779660, 88983512783, 345532572678, 1344372335524,
              5239988770268, 20457802016011 },
            // n = 3: https://oeis.org/A001931
            { 1, 1, 3, 15, 86, 534, 3481, 23502, 162913, 1152870, 8294738,
              60494549, 446205905, 3322769321, 24946773111, 188625900446,
              1435074454755, 10977812452428, 84384157287999, 651459315795897,
              5049008190434659, 39269513463794006, 306405169166373418 },
            // n = 4: https://oeis.org/A151830
            { 1, 1, 4, 28, 234, 2162, 21272, 218740, 2323730, 25314097,
              281345096, 3178474308, 36400646766, 421693622520, 4933625049464,
              58216226287844, 692095652493483 }
        }
    };
    return ((size_t) dim < cube_counts.size ()
            && (size_t) n_cells < cube_counts[dim].size ())
        ? cube_counts[dim][n_cells]
        : -1;
}

class Dim;
class Cells;
class PolyCube;
inline std::ostream& operator << (std::ostream&, Dim);
inline std::ostream& operator << (std::ostream&, const Cells&);
inline std::ostream& operator << (std::ostream&, const PolyCube&);

struct DimIterator;

#define CHECK_SIZE 0

struct Dim
{
    using value_t = int8_t;
    using vector_t = value_t __attribute__((vector_size(8)));
    using int_t = uint64_t;
    vector_t v = (vector_t) (int_t) 0;
    static inline constexpr int max_size = sizeof (vector_t) - 1;
    static inline constexpr vector_t all0 = (vector_t) (int_t) 0;
    Dim (std::initializer_list<value_t> &&il)
    {
        Dim::check_size ((int) il.size ());
        v = Dim::all0;
        set_size ((int) il.size ());
        int j = 0;
        for (auto i : il)
            v[j++] = i;
    }
    Dim (const vector_t &v) : v(v) {}

#if CHECK_SIZE
    static void check_size (int s) { assert (s >= 0 && s <= Dim::max_size); }
    void check_size (Dim d) const { assert (size () == d.size ()); }
#else
    static void check_size (int) {}
    void check_size (Dim) const {}
#endif

    static Dim zeros (int n_zeros)
    {
        Dim::check_size (n_zeros);
        Dim d{};
        d.set_size (n_zeros);
        return d;
    }
    int size () const
    {
        return v[Dim::max_size];
    }
    void set_size (int s)
    {
        Dim::check_size (s);
        v[Dim::max_size] = s;
    }
    int operator [] (int i) const
    {
        return v[i];
    }
    // Shift-invariant comparison, so (int_t) <=> (int_t) won't work.
    int cmp (Dim d) const
    {
        if (size () != d.size ())
            return size () - d.size ();
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
    void operator += (Dim d)
    {
        check_size (d);
        v += d.v;
        set_size (d.size ());
    }
    void operator -= (Dim d)
    {
        check_size (d);
        v -= d.v;
        set_size (d.size ());
    }
    Dim operator + (Dim d) const { Dim s (*this); s += d; return s; }
    Dim operator - (Dim d) const { Dim s (*this); s -= d; return s; }

    DimIterator begin () const;
    DimIterator end () const;
};

struct DimIterator
{
    Dim d;
    using iterator_category = std::forward_iterator_tag;
    DimIterator (Dim d) : d(d) {}
    bool operator != (DimIterator it) const
    {
        return d != it.d;
    }
    void operator ++ ()
    {
        for (int i = 0; ; ++i)
            if (d.v[i] == 1)
            {
                d.v[i] = -1;
                return;
            }
            else if (d.v[i] == -1)
            {
                d.v[i++] = 0;
                if (i < d.size ())
                    d.v[i] = 1;
                else
                    d.v = Dim::all0;
                return;
            }
    }
    Dim operator * () const
    {
        return d;
    }
};

inline DimIterator Dim::end () const { return DimIterator (Dim{}); }
inline DimIterator Dim::begin () const
{
    Dim d{};
    if (size ())
    {
        d.set_size (size ());
        d.v[0] = 1;
    }
    return DimIterator (d);
}
# if 0
struct Dim
{
    using value_type = int8_t;
    using vector_type = std::vector<value_type>;
    vector_type x;

    // ? Must not be unordered_set due to rehash ?
    using Pool = std::set<Dim>;
    static inline Pool pool;
    static void extend_pool (int dim, int n_cells);
    static const Dim* get (const Dim &d, bool may_update = 0)
    {
        const auto &f = pool.find (d);
        if (f != pool.end ())
            return & *f;
        else if (may_update)
        {
            assert (omp_get_thread_num () == 0);
            pool.insert (d);
            return &* pool.find (d);
        }
        else // ! may_update
        {
            std::cout << "Dim(" << d << ") not found\n";
            exit (2);
        }
    }

    Dim (const vector_type &x) : x(x) {}
    Dim (vector_type &&x) : x(std::move (x)) {}
    static Dim dim (int n)
    {
        return vector_type (n, 0);
    }
    int size () const
    {
        return (int) x.size ();
    }
    int cmp (const vector_type &y) const
    {
        if (x.size () != y.size ())
            return (int) x.size () - (int) y.size ();
        for (int i = 0; i < size (); ++i)
            if (x[i] != y[i])
                return x[i] - y[i];
        return 0;
    }
    int cmp (const Dim &d) const
    {
        return cmp (d.x);
    }
    int cmp (const Dim *d) const
    {
        return cmp (*d);
    }
    bool operator ==  (const Dim &d) const
    {
        return cmp (d) == 0;
    }
    bool operator <  (const Dim &d) const
    {
        return cmp (d) < 0;
    }
    ///////////////////////////////////////////////////////////////////////////
public:
    DimIterator begin (int reset_val = 0) const;
    DimIterator end () const;
    ///////////////////////////////////////////////////////////////////////////
};

struct DimIterator
{
    friend Dim;
    using iterator_category = std::forward_iterator_tag;
    using value_type        = Dim;
    using pointer           = value_type*;
    using reference         = value_type&;
    const Dim d0 {{}};
    Dim d {{}};
    bool done = false;
    int reset_value = 0;
    DimIterator (const Dim &d0, int reset_value = 0)
        : d0(d0), reset_value(reset_value)
    {
        d.x.resize (d0.x.size (), reset_value);
    }
private:
    DimIterator (bool done) : done(done) {}
public:
    DimIterator& operator ++ ()
    {
        for (size_t i = 0; i < d.x.size (); ++i)
            if (++d.x[i] >= d0.x[i] - reset_value)
                d.x[i] = (int8_t) reset_value;
            else
                return *this;
        done = true;
        return *this;
    }
    const Dim& operator * ()
    {
        return d;
    }

    bool operator == (const DimIterator &that) const
      {
          if (done || that.done)
              return done == that.done;
          return d == that.d;
      }
    bool operator != (const DimIterator &that) const
    {
        return !(*this == that);
    }
};

DimIterator Dim::begin (int reset_val) const
{
    return DimIterator (*this, reset_val);
}

DimIterator Dim::end () const
{
    return DimIterator (true);
}


inline void Dim::extend_pool (int dim, int n_cells)
{
    Dim d (Dim::vector_type (dim, n_cells));
    for (auto &&it = d.begin(-2); !it.done; ++it)
    {
        bool neg_p = false;
        int dist = 0;
        for (int i : (*it).x)
        {
            neg_p |= i < 0;
            dist += i;
        }
        if (neg_p || (! neg_p && dist <= n_cells))
            (void) Dim::get (*it, true);
    }
}
#endif

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
                if (d > *p)
                    continue;
                if (d < *p)
                    cells.emplace (p, d);
                break;
            }
    }
    void erase (Dim d)
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
    int cmp (const Cells &c) const
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
    unsigned hash () const
    {
        unsigned h = 0;
        for (Dim d : cells)
            h = h * 13 + (Dim::int_t) d.v;
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
            if (i >= c.size ())
            {
                std::cout << "change " << i << " of " << c << "\n";
                exit (1);
            }
            assert (i < c.size ());
            c.v[i] += off;
        }
    }
};


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
    using Vector = std::vector<MuxSet>; // Indexed by corona size.
    unsigned m_hash = 0;
    Cells m_cubes;

    const Cells corona () const
    {
        Cells cora;
        for (Dim d : m_cubes.cells)
            for (Dim delta : d)
                if (! m_cubes.contains (d + delta))
                    cora.add (d + delta);
        return cora;
    }
    int corona_size () const
    {
        return corona ().size ();
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
    unsigned hash () const
    {
        return m_hash;
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

    // Way 3
    void add_sprouts_way3 (Vector &vms) const
    {
        for (Dim d : corona().cells)
        {
            PolyCube pc (*this);
            pc.add (d);
            auto &slot = vms[pc.corona_size ()];
            slot.mux.lock ();
            slot.set.emplace (std::move (pc));
            slot.mux.unlock ();
        }
    }

    // Way 3
    static void add_sprouts_way3 (int dim, int n,
                                  Vector &vset2, const Vector &vset)
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
#pragma omp parallel for schedule(runtime)
        for (size_t j = 0; j < vpc.size (); ++j)
            vpc[j]->add_sprouts_way3 (vset2);
    }

    // Way 4
    void add_sprouts_way4 (Vector &vms) const
    {
        for (Dim d : corona().cells)
        {
            PolyCube pc (*this);
            pc.add (d);
            MuxSet &slot = vms[pc.hash () % vms.size ()];
            slot.mux.lock ();
            slot.set.emplace (std::move (pc));
            slot.mux.unlock ();
        }
    }

    // Way 4
    static void add_sprouts_way4 (int n_slots,
                                  Vector &vset2, const Vector &vset)
    {
        Vector v (n_slots);
        vset2.swap (v); // Since resize() doesn't like std::mutex.

#pragma omp parallel for schedule(runtime)
        for (size_t j = 0; j < vset.size (); ++j)
        {
            for (const auto &pc : vset[j].set)
                pc.add_sprouts_way4 (vset2);
        }
    }

    // Way 6, 7
    void add_sprouts_merge (Set &set) const
    {
        for (Dim d : corona().cells)
        {
            PolyCube pc (*this);
            pc.add (d);
            PolyCube::merge (set, std::move (pc));
        }
    }

    // Way 6
    static void add_sprouts_6reduce (Set &set2, const Set &set)
    {
        std::vector<const PolyCube*> pc;
        pc.resize (set.size ());
        int j = 0;
        for (const auto &p : set)
            pc[j++] = &p;

#pragma omp parallel for schedule(runtime)
        for (size_t j = 0; j < pc.size (); ++j)
        {
            Set s;
            pc[j]->add_sprouts_merge (s);
#pragma omp critical
            PolyCube::merge (set2, s);
        }
    }

    // Way 7
    static void add_sprouts_7reduce (Set &set2, const Set &set, int n_pc)
    {
        std::vector<const PolyCube*> pc;
        pc.resize (set.size ());
        int j = 0;
        for (const auto &p : set)
            pc[j++] = &p;

        int dim = set.begin() -> m_cubes.cells.front().size ();
        int n_cells = 1 + (int) set.begin () -> m_cubes.size ();
        const int64_t n_cubes = cube_count (dim, n_cells);
#pragma omp parallel for schedule(runtime)
        for (size_t j = 0; j < pc.size (); j += n_pc)
        {
            Set s;
            bool done = false;
            for (int k = 0; k < n_pc && !done; ++k)
                if (j + k < pc.size ())
                    pc[j + k]->add_sprouts_merge (s);
#pragma omp critical
            {
                PolyCube::merge (set2, s);
                if (n_cubes > 1 && n_cubes == (int64_t) set2.size ())
                    done = true;
            }
        }
    }

    // Univariate polynomial over Z in sparse representation.
    using Poly = std::map<int,int>;

    static Poly get_poly (const Set &set)
    {
        Poly poly;
        for (const auto &pc : set)
        {
            const int coro = pc.corona_size ();
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

    static void add_poly (Poly &p, const Poly &q)
    {
        for (const auto &q_mono : q)
        {
            const int expo = q_mono.first;
            const int coef = q_mono.second;
            const auto &p_mono = p.find (expo);
            if (p_mono == p.end ())
                p[expo] = coef;
            else
                p_mono->second += coef;
        }
    }

    static Poly get_poly_way4 (const Vector &vms)
    {
        Poly poly;
#pragma omp parallel for
        for (size_t j = 0; j < vms.size (); ++j)
        {
            Poly p = PolyCube::get_poly (vms[j].set);
#pragma omp critical
            PolyCube::add_poly (poly, p);
        }
        return poly;
    }

    enum { POLY_TEX, POLY_LIST };

    static void print_poly (int n, const Poly &poly,
                            int style, const char *var = "q")
    {
        if (style == POLY_TEX)
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
        }
        else if (style == POLY_LIST)
        {
            bool start = true;
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
    }
};

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

inline std::ostream& operator << (std::ostream &ost, const Cells &c)
{
    ost << "{#" << c.size ();
    for (auto c : c.cells)
        ost << " " << c;
    return ost << " }";
}

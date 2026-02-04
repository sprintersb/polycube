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
#include <atomic>

#include <cstdint>
#include <cinttypes>
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
#error DIM=?
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
    Dim operator + (Dim d) const { Dim s (*this); s += d; return s; }
    Dim operator - (Dim d) const { Dim s (*this); s -= d; return s; }
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

    DimIterator begin () const;
    DimIterator end () const;
};

struct Box
{
    Dim lo, hi;
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
        cells.resize (1 + cells.size (), Dim{});

        for (size_t i = 0; ; ++i)
            if (i == cells.size () - 1)
            {
                cells[i] = d;
                break;
            }
            else
            {
                const int c = d.cmp (cells[i]);
                assert (c != 0 && "Assume we always increase cells");
                if (c > 0)
                    continue;
                if (c < 0)
                {
                    for (size_t j = cells.size () - 1; j > i; --j)
                        cells[j] = cells[j - 1];
                    cells[i] = d;
                }
                break;
            }
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
            if (i >= c.size ())
            {
                std::cout << "change " << i << " of " << c << "\n";
                exit (1);
            }
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
};


struct Corona
{
    struct DimHash
    {
        hash_t operator () (Dim d) const
        {
            return d.hash ();
        }
    };
    std::unordered_set<Dim, DimHash> cells;

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
    using Vector = std::vector<MuxSet>; // Indexed by corona size.
    hash_t m_hash = 0;
    Cubes m_cubes;

    const Corona corona () const
    {
        Corona cora;
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
    hash_t hash () const
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

    // Way 4
    int add_sprouts_way4 (Vector &vms) const
    {
        int new_count = 0;
        for (Dim d : corona().cells)
        {
            PolyCube pc (*this);
            pc.add (d);
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
                                  Vector &vset2, const Vector &vset)
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
                pc_count += pc.add_sprouts_way4 (vset2);

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
                const int coro = pc.corona_size ();
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

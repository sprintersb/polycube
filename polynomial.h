// -*- c++ -*-

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
        (void) x;
        Assert (x == 1, "todo: evaluate Poly at x != 1");
        int64_t y = 0;
        for (const auto &mono : a_)
            y += mono.second;
        return y;
    }

#pragma omp declare                             \
    reduction (+ : Poly : omp_out += omp_in)    \
    initializer (omp_priv = omp_orig)

#pragma omp declare                                                 \
    reduction (vecadd : std::vector<T> : vecadd (omp_out, omp_in))  \
    initializer (omp_priv = omp_orig)

    static void vecadd (std::vector<T> &a, const std::vector<T> &b)
    {
        Assert (a.size () == b.size (), "vector add of different sizes");
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
                Assert (coro < (int) v.size (), "impossible corona size");
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
                Assert (coro < (int) v.size (), "impossible corona size");
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

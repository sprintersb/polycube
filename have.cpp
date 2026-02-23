#include "polycube.h"
#include "diagnostic.h"

int64_t PolyCube::Have::get_fixed_count () const
{
    int64_t cnt = 0;

    if (fixed.count)
        return * fixed.count;
    else if (fixed.has_set ())
        return fixed.set.size ();
    else if (fixed.has_vec ())
        for (const auto &ms : fixed.vec)
            cnt += ms.set.size ();
    else if (corona_polynomial)
        return (*corona_polynomial) (1);
    else if (Cubes::can_canonicalize && free.has_vec ())
#pragma omp parallel for schedule(dynamic) reduction(+ : cnt)
        for (size_t i = 0; i < free.vec.size (); ++i)
        {
            int64_t c = 0;
            for (const auto &pc : free.vec[i].set)
                c += pc.multiplicity ();
            cnt += c;
        }
    else if (Cubes::can_canonicalize && free.has_set ())
        for (const auto &pc : free.set)
            cnt += pc.multiplicity ();
    else
        warning ("don't know how to obtain fixed cubes count");

    return cnt ? cnt : -1;
}

int64_t PolyCube::Have::get_free_count () const
{
    int64_t cnt = 0;

    if (free.count)
        return * free.count;
    else if (free.has_set ())
        return free.set.size ();
    else if (fixed.has_vec ())
        for (const auto &ms : fixed.vec)
            cnt += ms.set.size ();
    else if (Cubes::can_canonicalize && fixed.has_vec ())
#pragma omp parallel for schedule(dynamic) reduction(+ : cnt)
        for (size_t i = 0; i < fixed.vec.size (); ++i)
        {
            int64_t c = 0;
            for (const auto &pc : fixed.vec[i].set)
                c += pc.is_canonical ();
            cnt += c;
        }
    else if (Cubes::can_canonicalize && fixed.has_set ())
        for (const auto &pc : fixed.set)
            cnt += pc.is_canonical ();
    else
        warning ("don't know how to obtain free cubes count");

    return cnt ? cnt : -1;
}

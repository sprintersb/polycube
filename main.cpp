#include <vector>
#include <iostream>
// C
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cinttypes>
#include <cmath>
// Own
#include "polycube.h"

int main_polycube (int argc, char *argv[])
{
    int dim = 2;
    int level = 10;
    int way = 0;
    int extra_spice = 0;
    int leap = 1;

    if (argc > 1)   sscanf (argv[1], "%i", &dim);
    if (argc > 2)   sscanf (argv[2], "%i", &level);
    if (argc > 3)   sscanf (argv[3], "%i", &way);
    if (argc > 4)   sscanf (argv[4], "%i", &extra_spice);
    if (argc > 5)   sscanf (argv[5], "%i", &leap);

    PolyCube::canonical_only = way > 0 && (2 <= dim && dim <= 4);

#if defined CUBES_ARRAY
    std::cout << "CUBES_ARRAY  DIM=" << DIM << "  CELLS=" << CELLS << "\n";
#else
    std::cout << "CUBES_VECT DIM=" << DIM << "\n";
#endif

    assert (way == 0 || way == 4 || way == 5);
    assert (dim == DIM);
#ifdef CUBES_ARRAY
    assert (level <= CELLS);
#endif
    if (way == 5)
        assert (level - leap >= 1);

    const int max_threads = omp_get_max_threads ();
    std::cout << "Max threads: " << max_threads << "\n";

    // See birthday paradox.
    const double p_slow_collide = 0.2;
    const int n_slots = (int) (max_threads * max_threads / 2 / p_slow_collide);
    std::cout << "Slots      : " << n_slots << "\n";

    max_possible_corona = 2 * (dim - 1) * level + 2;
    std::cout << "maxi corona: " << max_possible_corona << "\n";
    std::cout << "canonical o: " << PolyCube::canonical_only << "\n";

    std::vector<PolyCube::Set> set (1 + level);     // Way 0
    std::vector<PolyCube::Vector> vset (1 + level); // Way 4, 5
    std::vector<int> smallest_corona (1 + level, 0);

    for (int i = 1; i <= level; ++i)
    {
        std::cout << "== " << i << " ==\n";
        PolyCube::Poly poly;

        if (i == 1)
        {
            PolyCube pc1 (nullptr, Dim::all (0));
            if (way == 4 || way == 5)
            {
                // Index is hash % n_slots.
                PolyCube::Vector v (n_slots);
                vset[1].swap (v);
                vset[1][pc1.hash () % n_slots].set.emplace (pc1);
                if (way == 5)
                    poly = PolyCube::Poly (vset[1]);
            }
            else
                set[1].emplace (pc1);
        }
        else
        {
            stat[0] = 0;
            stat[1] = 0;
            stat[2] = 0;
            if (way == 4)
            {
                int corona_margin = extra_spice;
                int max_corona = corona_margin > 0
                    ? corona_margin + smallest_corona[i - 1]
                    : -1;
                PolyCube::Filter filter =
                    [max_corona](const PolyCube &pc)
                    {
                        return ! pc.has_large_corona (max_corona);
                    };
                if (max_corona <= 0)
                    filter = nullptr;
                PolyCube::add_sprouts_way4 (i, n_slots, 0, vset[i], vset[i - 1],
                                            200, filter, nullptr);
            }
            else if (way == 5)
            {
                if (i <= level - leap)
                    poly = PolyCube
                        ::get_sprouts_poly_way5 (i, n_slots,
                                                 extra_spice, 0,
                                                 vset[i], vset[i - 1]);
                else if (i == level)
                    poly = PolyCube
                        ::get_sprouts_poly_way5 (i, n_slots,
                                                 extra_spice, leap,
                                                 vset[i], vset[i - leap]);
                else
                {
                    std::cout << "leaped\n";
                    continue;
                }
            }
            else if (way == 0)
            {
                for (const auto &pc : set[i - 1])
                    if (way == 0) pc.add_sprouts (set[i]);
            }
        }

        // Stat about fraction of fast canonicalization.
        auto tot = stat[0] + stat[1];
        double f0 = tot ? (double) +stat[0] / tot : 0.0;
        double f1 = tot ? (double) +stat[1] / tot : 0.0;
        printf ("Stat: %" PRIi64 "(%.2f%%)  %" PRIi64 "(%.2f%%)\n",
                +stat[0], 100 * f0, +stat[1], 100 * f1);

        if (way == 4)
            poly = PolyCube::Poly (vset[i]);
        else if (way == 0)
            poly = PolyCube::Poly (set[i]);

        uint64_t ccount = poly (1);

        smallest_corona[i] = poly.a_.begin ()->first;
        std::cout << ccount << " polycubes"
                  << "  (coro min: " << smallest_corona[i] << ")\n";
        if (way == 4 && extra_spice)
        {
            const double r = 2 * sqrt (2 * i - 1) + 2;
            printf ("coro min / calculated = %d / %.2f = %.2f\n",
                    smallest_corona[i], r, smallest_corona[i] / r);
        }

        poly.print (i, PolyCube::Poly::POLY_TEX);
        std::cout.flush();

        if (way == 4 && extra_spice > 0)
        {
            PolyCube::Set small_corona = PolyCube::find_min_corona (vset[i]);
            if (! small_corona.empty ())
            {
                const auto &best = * small_corona.begin();
                std::cout << best.ascii ();
                std::cout << BorderFinder (best).border().svg() << "\n";
            }
            std::cout.flush();
        }

        if (way != 4 || extra_spice /* corona_margin */ <= 0)
            if (cube_count (dim, i) >= 0)
                assert ((int64_t) ccount == cube_count (dim, i)
                        && "verify polycube count");

        if (0 && i <= 3 && DIM == 2)
            for (const auto &ms: vset[i])
                for (const auto &pc : ms.set)
                {
                    auto &&ps = BorderFinder (pc).border();
                    for (const auto &pgon : ps)
                    {
                        std::cout << pgon;
                        std::cout << pc.ascii ();
                    }
                }
    }
    exit (0); // Faster than waiting for all them destructors.
}


int main (int argc, char *argv[])
{
    return main_polycube (argc, argv);
}

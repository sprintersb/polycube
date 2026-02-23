#include <vector>
#include <iostream>
#include <optional>
// C
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cinttypes>
#include <cmath>
// Own
#include "cubes-border.h"
#include "polycube.h"
#include "diagnostic.h"
#include "is-a-tty.h"

[[gnu::constructor]]
static void sync_stdio ()
{
    // Allow mixing std::cout and stdout and keep their order.
    std::ios::sync_with_stdio (true);
}

void show_growth_rates ()
{
    printf ("Fixed cubes growth rates\n");
    for (int dim = 2; dim <= 5; ++dim)
    {
        printf ("dim = %d:", dim);
        for (int i = 5; ; ++i)
            if (auto a = cube_count (dim, i); a > 0)
                if (auto b = cube_count (dim, i - 1); b > 0)
                    printf ("  %.2f", (double) a / b);
                else break;
            else break;
        printf ("\n");
    }

    printf ("Free cubes growth rates\n");
    for (int dim = 2; dim <= 5; ++dim)
    {
        printf ("dim = %d:", dim);
        for (int i = 5; ; ++i)
            if (auto a = cube_count_canon (dim, i); a > 0)
                if (auto b = cube_count_canon (dim, i - 1); b > 0)
                    printf ("  %.2f", (double) a / b);
                else break;
            else break;
        printf ("\n");
    }

    printf ("Fixed cubes / free cubes ratio\n");
    for (int dim = 2; dim <= 5; ++dim)
    {
        printf ("dim = %d:", dim);
        for (int i = 5; ; ++i)
            if (auto a = cube_count (dim, i); a > 0)
                if (auto b = cube_count_canon (dim, i); b > 0)
                    printf ("  %.2f", (double) a / b);
                else break;
            else break;
        printf ("\n");
    }
}

void print_stat ()
{
    if (!Cubes::take_stat)
        return;

    // Stat about fraction of fast canonicalization.
    const int64_t s0 = stat[STAT_FAIL];
    const int64_t s1 = stat[STAT_SUCC];
    const int64_t s2 = stat[STAT_COST_FAIL];
    const int64_t s3 = stat[STAT_COST_SUCC];
    auto tot = s0 + s1;
    if (tot)
    {
        const double f0 = tot ? (double) s0 / tot : 0.0;
        const double f1 = tot ? (double) s1 / tot : 0.0;
        printf ("Stat: %" PRIi64 "(%.2f%%)  %" PRIi64 "(%.2f%%)",
                s0, 100 * f0, s1, 100 * f1);
        const double cf0 = s0 ? (double) s2 / s0 : 0.0;
        const double cf1 = s1 ? (double) s3 / s1 : 0.0;
        printf ("\t Cost factor: %.1f + %.1f = %.1f\n",
                cf0, cf1, f0 * cf0 + f1 * cf1);
    }
    else
        printf ("Cost factor: %.1f\n", 0. + gjl::hyperoctahedral_order(DIM));
}


int main_polycube (int argc, char *argv[])
{
    show_growth_rates ();

    int dim = 2;
    int level = 10;
    int way = 0;
    int extra_spice = 0;
    int leap = 1;
    Cubes::take_stat = false;

    if (argc > 1)   sscanf (argv[1], "%i", &dim);
    if (argc > 2)   sscanf (argv[2], "%i", &level);
    if (argc > 3)   sscanf (argv[3], "%i", &way);
    if (argc > 4)   sscanf (argv[4], "%i", &extra_spice);
    if (argc > 5)   sscanf (argv[5], "%i", &leap);

    PolyCube::canonical_only = Cubes::can_canonicalize;

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

    const int max_threads = omp_get_max_threads ();
    std::cout << "Max threads: " << max_threads << "\n";

    // See birthday paradox.
    const double p_slot_collide = 0.2;
    const int n_slots
        = std::max<int> (1, max_threads * max_threads / 2 / p_slot_collide);
    std::cout << "Slots      : " << n_slots << "\n";

    max_possible_corona = 2 * (dim - 1) * level + 2;
    std::cout << "maxi corona: " << max_possible_corona << "\n";
    std::cout << "can canonic: " << Cubes::can_canonicalize << "\n";
    std::cout << "take stat  : " << Cubes::take_stat << "\n";

    std::vector<PolyCube::Have> have (1 + level);

    auto &want = PolyCube::want;
    want.free.cubes = Cubes::can_canonicalize;
    want.free.count = Cubes::can_canonicalize;
    want.fixed.cubes = ! Cubes::can_canonicalize;
    want.fixed.count = true;
    want.n_parts = extra_spice;

    want.progress = is_a_tty (stdout);
    Progress<int64_t>::quiet = ! want.progress;
    int progress_at = 200;

    want.corona_polynomial = 0;//true;

    if (want.n_parts < 2)
        leap = 0;
    leap = std::min (leap, level - 1);

    for (int i = 1; i <= level; ++i)
    {
        std::cout << "== " << i << " ==\n";

        want.n_cells = i;
        want.n_cells_final = level;

        auto &poly = have[i].corona_polynomial;

        if (i == 1)
        {
            PolyCube pc1 (nullptr, Dim::all (0));
            if (way == 4 || way == 5)
            {
                // Index is hash % n_slots.
                PolyCube::Vector v (n_slots);
                have[i]().vec.swap (v);
                have[i]().vec[pc1.hash () % n_slots].set.emplace (pc1);
                if (way == 5)
                    poly = PolyCube::Poly (have[i]().vec);
            }
            else
                have[i]().set.emplace (pc1);
            have[i]().count = 1;
            have[i].smallest_corona = pc1;
        }
        else
        {
            for (auto &s : stat)
                s = 0;

            if (way == 4)
            {
                int corona_margin = extra_spice;
                int max_corona = corona_margin > 0
                    ? (corona_margin
                       + have[i - 1].smallest_corona->corona().size())
                    : -1;
                PolyCube::Filter filter = nullptr;
                if (max_corona > 0)
                    filter = [max_corona](const PolyCube &pc)
                    {
                        return ! pc.has_large_corona (max_corona);
                    };
                have[i]().count = PolyCube::sprout_way4 (n_slots, progress_at,
                                                         have[i], have[i - 1],
                                                         filter, nullptr);
            }
            else if (way == 5)
            {
                want.leap = i == level && leap ? leap : 0;
                if (want.leap)
                    want.free.cubes = want.fixed.cubes = false;
                const auto &hin = have[want.leap ? i - want.leap : i - 1];
                if (i <= level - leap || i == level)
                    PolyCube::sprout_way5 (n_slots, progress_at, have[i], hin);
                else
                {
                    std::cout << "leaped\n";
                    continue;
                }
            }
            else if (way == 0)
            {
                for (const auto &pc : have[i - 1]().set)
                    pc.add_sprouts (have[i]().set);
                have[i]().count = have[i]().set.size ();
            }
        }

        print_stat ();

        if (int ci = cube_count_canon (dim, i); ci >= 0)
            if (int ci1 = cube_count_canon (dim - 1, i); ci1 >= 0)
                printf ("free %dd = %.1f%%\n", DIM - 1, 100. * ci1 / ci);

        if (want.corona_polynomial)
        {
            if (way == 4)
                poly = PolyCube::Poly (have[i]().vec);
            else if (way == 0)
                poly = PolyCube::Poly (have[i]().set);
        }

        int64_t ccount = have[i].get_fixed_count ();

        std::cout << ccount << " polycubes";
        if (have[i].smallest_corona)
            std::cout << "  (coro min: "
                      << have[i].smallest_corona->corona().size() << ")";
        std::cout << "\n";
        if (way == 4 && extra_spice && have[i].smallest_corona)
        {
            const int scs = have[i].smallest_corona->corona().size();
            const double r = 2 * sqrt (2 * i - 1) + 2;
            printf ("coro min / calculated = %d / %.2f = %.2f\n",
                    scs, r, scs / r);
        }

        if (poly)
        {
            poly->print (i, PolyCube::Poly::POLY_TEX);
            std::cout.flush();
        }

        if (way == 4 && extra_spice > 0)
        {
            auto small_corona = PolyCube::find_min_corona (have[i]().vec);
            if (! small_corona.empty ())
            {
                const auto &best = * small_corona.begin();
                std::cout << best.ascii ();
                std::cout << cubes_border (best).svg() << "\n";
            }
            std::cout.flush();
        }

        if (way != 4 || extra_spice /* corona_margin */ <= 0)
            if (cube_count (dim, i) >= 0
                && ccount != cube_count (dim, i))
            {
                error ("cube count %" PRIi64 " != %" PRIi64 " expected count",
                       ccount, cube_count (dim, i));
            }

        if (0 && i <= 3 && DIM == 2)
            for (const auto &ms: have[i]().vec)
                for (const auto &pc : ms.set)
                {
                    auto &&ps = cubes_border (pc);
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

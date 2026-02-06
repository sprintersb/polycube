#include "polycube.h"
#include <cstdlib>

int main_polycube (int argc, char *argv[])
{
    int dim = 2;
    int level = 10;
    int way = 0;
    int n_pc = 1;

    if (argc > 1)   sscanf (argv[1], "%i", &dim);
    if (argc > 2)   sscanf (argv[2], "%i", &level);
    if (argc > 3)   sscanf (argv[3], "%i", &way);
    if (argc > 4)   sscanf (argv[4], "%i", &n_pc);

    assert (way == 0 || way == 4);
    assert (way != 5 && way != 1 && way != 2);

    std::vector<PolyCube::Set> set (1 + level);     // Way 0
    std::vector<PolyCube::Vector> vset (1 + level); // Way 4

    for (int i = 1; i <= level; ++i)
    {
        std::cout << "== " << i << " ==\n";

        if (i == 1)
        {
            PolyCube pc1 (nullptr, Dim::zeros (dim));
            if (way == 4)
            {
                // Index is hash % n_pc.
                PolyCube::Vector v (n_pc);
                vset[1].swap (v);
                vset[1][pc1.hash () % n_pc].set.emplace (pc1);
            }
            else
                set[1].emplace (pc1);
        }
        else
        {
            if (way == 4)
                PolyCube::add_sprouts_way4 (n_pc, vset[i], vset[i - 1]);
            else if (way == 0)
            {
                for (const auto &pc : set[i - 1])
                    if (way == 0) pc.add_sprouts (set[i]);
            }
        }

        uint64_t ccount = -1;
        PolyCube::Poly poly;

        if (way == 3 || way == 4)
        {
            int n_polycubes = 0;
            for (const auto &ms : vset[i])
                n_polycubes += ms.set.size ();
            std::cout << (ccount = n_polycubes) << " polycubes\n";
            if (way == 3)
                poly = PolyCube::get_poly (vset[i]);
            else if (way == 4)
                poly = PolyCube::get_poly_way4 (vset[i]);
        }
        else
        {
            std::cout << (ccount = set[i].size()) << " polycubes\n";
            poly = PolyCube::get_poly (set[i]);
        }

        PolyCube::print_poly (i, poly, PolyCube::POLY_TEX);

        if (cube_count (dim, i) >= 0)
            assert ((int64_t) ccount == cube_count (dim, i)
                    && "verify polycube count");
    }
    exit (0); // Faster than waiting for all them destructors.
}


int main (int argc, char *argv[])
{
    return main_polycube (argc, argv);
}

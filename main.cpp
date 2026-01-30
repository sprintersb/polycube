#include "polycube.h"

int main_polycube (int argc, char *argv[])
{
    int dim = 2;
    int level = 10;
    int way = 0;

    if (argc > 1)   sscanf (argv[1], "%i", &dim);
    if (argc > 2)   sscanf (argv[2], "%i", &level);
    if (argc > 3)   sscanf (argv[3], "%i", &way);

    std::vector<PolyCube::Set> set (1 + level);     // Way 0..2
    std::vector<PolyCube::Vector> vset (1 + level); // Way 3

    for (int i = 1; i <= level; ++i)
    {
        std::cout << "== " << i << " ==\n";

        if (i == 1)
        {
            PolyCube pc1;
            pc1.add (Dim::dim (dim));
            if (way == 3)
            {
                PolyCube::Vector v (1 + 2 * dim);
                vset[1].swap (v);
                assert (pc1.corona.size () == 2 * dim);
                vset[1][2 * dim].set.emplace (pc1);
            }
            else
                set[1].emplace (pc1);
        }
        else
        {
            if (way == 2)
                PolyCube::add_sprouts (set[i], set[i - 1]);
            else if (way == 3)
            {
                PolyCube::add_sprouts (dim, i, vset[i], vset[i - 1]);
            }
            else
                for (const auto &pc : set[i - 1])
                {
                    if (way == 0) pc.add_sprouts (set[i]);
                    if (way == 1) pc.add_sprouts_omp (set[i]);
                }
        }

        if (way == 3)
        {
            int n_polycubes = 0;
            for (const auto &ms : vset[i])
                n_polycubes += ms.set.size ();
            std::cout << n_polycubes << " polycubes\n";
            auto poly = PolyCube::get_poly (vset[i]);
            PolyCube::print_poly (i, poly);
        }
        else
        {
            std::cout << set[i].size() << " polycubes\n";
            auto poly = PolyCube::get_poly (set[i]);
            PolyCube::print_poly (i, poly);
        }
    }
    return 0;
}


int main (int argc, char *argv[])
{
    return main_polycube (argc, argv);
}

#include "polycube.h"

int main_polycube (int argc, char *argv[])
{
    int dim = 2;
    int level = 10;
    int way = 0;

    if (argc > 1)   sscanf (argv[1], "%i", &dim);
    if (argc > 2)   sscanf (argv[2], "%i", &level);
    if (argc > 3)   sscanf (argv[3], "%i", &way);

    std::vector<PolyCube::Set> set;
    set.push_back (PolyCube::Set {});

    for (int i = 1; i <= level; ++i)
    {
        std::cout << "== " << i << " ==\n";
        set.push_back (PolyCube::Set {});
        PolyCube::Set &s = set[i];
        if (i == 1)
        {
            PolyCube pc1;
            pc1.add (Dim::dim (dim));
            s.emplace (pc1);
        }
        else
            for (const auto &pc : set[i - 1])
            {
                if (way == 0) pc.add_sprouts (s);
                if (way == 1) pc.add_sprouts_omp (s);
            }
        std::cout << s.size() << " polycubes\n";
        auto poly = PolyCube::get_poly (s);
        PolyCube::print_poly (i, poly);
    }
    return 0;
}


int main (int argc, char *argv[])
{
    return main_polycube (argc, argv);
}

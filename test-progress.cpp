#include <unistd.h>
#include "progress.h"
#include <iostream>
#include <cstdint>

using Pro = Progress<int64_t>;

Pro::Printer pp = [] (std::ostringstream &ost, Pro &p, int i) -> void
{
    ost << i << " Cubs = " << 100.0 * i / p.total << "%";
};

Pro::Printer pbar = [] (std::ostringstream &ost, Pro &p, int i) -> void
{
    ost << i << ": ";
    p.print_bar (ost, 79, (double) i / p.total);
};

void prog ()
{
    std::cout << "Test progress\n";
    Pro p (1000, 1, pp);
    for (int i = 0; i < 1000; ++i)
    {
        usleep (2000);
        p.update (i);
    }
    p.done ();

    std::cout << "Test progress bar\n";
    p.printer = & pbar;
    p.reset ();
    for (int i = 0; i < 1000; ++i)
    {
        usleep (2000);
        p.update (i);
    }
}

int main ()
{
    prog();
}

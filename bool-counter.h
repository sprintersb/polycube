// -*- c++ -*-
#ifndef BOOL_COUNTER_H
#define BOOL_COUNTER_H

#include <iostream>
#include <cstdio>

struct BoolCounter
{
    int c[2] = { 0, 0 };
    void reset ()
    {
        c[0] = c[1] = 0;
    }
    bool operator += (bool b)
    {
        ++c[b];
        return b;
    }
    friend std::ostream& operator << (std::ostream &ost, const BoolCounter &c);
};

inline std::ostream& operator << (std::ostream &ost, const BoolCounter &c)
{
    ost << c.c[0] << "/" << c.c[1];
    if (c.c[1])
    {
        char part[20];
        sprintf (part, "%.2f%%", 100.0 * c.c[0] / c.c[1]);
        ost << " = " << part;
    }
    return ost;
}

#endif // BOOL_COUNTER_H

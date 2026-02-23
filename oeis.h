// -*- c++ -*-
#ifndef OEIS_H
#define OEIS_H

#include <array>
#include <vector>
#include <cstdint>

namespace oeis
{
    using Vector = std::vector<int64_t>;

    struct Sequence : Vector
    {
        const char *const name;
        using T = Vector::value_type;
        using value_type = T;
        Sequence () : Sequence ("Void", {}) {}
        Sequence (const char *name, std::vector<T> v) : Vector(v), name(name)
        {}
        T at (int i) const
        {
            return i < (int) size () ? Vector::at (i) : -1;
        }
        T operator [] (int i) const
        {
            return at (i);
        }
    };

    inline const Sequence empty;

    inline const std::array<Sequence, 5 + 1> seq_fixed =
    {
        empty,
        empty,
        { "A001168", // n = 2: https://oeis.org/A001168
          { 1, 1, 2, 6, 19, 63, 216, 760, 2725, 9910, 36446, 135268, 505861,
            1903890, 7204874, 27394666, 104592937, 400795844, 1540820542,
            5940738676, 22964779660, 88983512783, 345532572678, 1344372335524,
            5239988770268, 20457802016011 } },
        { "A001931", // n = 3: https://oeis.org/A001931
          { 1, 1, 3, 15, 86, 534, 3481, 23502, 162913, 1152870, 8294738,
            60494549, 446205905, 3322769321, 24946773111, 188625900446,
            1435074454755, 10977812452428, 84384157287999, 651459315795897,
            5049008190434659, 39269513463794006, 306405169166373418 } },
        { "A151830", // n = 4: https://oeis.org/A151830
          { 1, 1, 4, 28, 234, 2162, 21272, 218740, 2323730, 25314097,
            281345096, 3178474308, 36400646766, 421693622520, 4933625049464,
            58216226287844, 692095652493483 } },
        { "A151831", // n = 5: https://oeis.org/A151831
          { 1, 1, 5, 45, 495, 6095, 80617, 1121075, 16177405, 240196280,
            3648115531, 56440473990, 886696345225, 14111836458890,
            227093585071305, 3689707621144614 } },
    };

    inline const std::array<Sequence, 5 + 1> seq_free =
    {
        empty,
        empty,
        { "A000105", // n = 2: https://oeis.org/A000105
          { 1, 1, 1, 2, 5, 12, 35, 108, 369, 1285, 4655, 17073, 63600,
            238591, 901971, 3426576, 13079255, 50107909, 192622052,
            742624232, 2870671950, 11123060678, 43191857688, 168047007728,
            654999700403, 2557227044764, 9999088822075, 39153010938487,
            153511100594603 } },
        { "A038119", // n = 3: https://oeis.org/A038119
          { 1, 1, 1, 2, 7, 23, 112, 607, 3811, 25413, 178083, 1279537,
            9371094, 69513546, 520878101, 3934285874, 29915913663,
            228779330204, 1758309223457, 13573319825615, 105192814197984,
            818136047201932, 6383528588447574 } },
        { "A068870", // n = 4: https://oeis.org/A068870
          { 1, 1, 1, 2, 7, 26, 147, 1019, 8699, 82535, 846042, 9078720,
            /* GJL */ 100651853, 1141767844,
            /* GJL */ 13177518932 /* 26000 min = 18 days, 3 parts */ } },
        { "", // n = 5: ???
          { /* GJL */ 1, 1, 1, 2, 7, 26, 153, 1123, 10708, 119120,
            /*`GJL */ 1493722, 20252600, 290460057,
            /* GJL */ 4335535057 /* 44000 min = 31 days, 2 parts */ } },
    };

    inline const Sequence& cubes_fixed (int dim)
    {
        return dim < (int) seq_fixed.size () ? seq_fixed[dim] : empty;
    }
    inline Sequence::value_type cubes_fixed (int dim, int cells)
    {
        return cubes_fixed (dim)[cells];
    }

    inline const Sequence& cubes_free (int dim)
    {
        return dim < (int) seq_free.size () ? seq_free[dim] : empty;
    }
    inline Sequence::value_type cubes_free (int dim, int cells)
    {
        return cubes_free (dim)[cells];
    }
} // oeis
#endif // OEIS_H

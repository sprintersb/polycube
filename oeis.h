// -*- c++ -*-
#ifndef OEIS_H
#define OEIS_H

#include <vector>
// C
#include <cstdint>
#include <cstring>
// Own
#include "diagnostic.h"

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
    using Sequences = std::vector<Sequence>;

    inline const Sequence empty;

    inline const Sequences seq_fixed =
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

    inline const Sequences seq_free =
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

    inline const Sequences seq_symm =
    {
        empty,
        empty,
        { "A030227", // n = 2: https://oeis.org/A030227
          { 1, 1, 1, 2, 3, 6, 10, 20, 34, 70, 121, 250, 441, 912, 1630, 3375,
            6092, 12624, 22961, 47616, 87136, 180811, 332549, 690398, 1275166,
            2648422, 4909364, 10199792, 18966700, 39416488, 73497642,
            152777230, 285569898, 593717419, 1112188817, 2312672439,
            4340728280, 9027238683, 16973536668, 35303017659 } },
        { "A007743", // n = 3: https://oeis.org/A007743
          { 1, 1, 1, 2, 6, 17, 58, 191, 700, 2515, 9623, 36552, 143761, 564443,
            2259905, 9057278, 36705846, 149046429, 609246350, 2495727647,
            10267016450, 42322763940, 174974139365 } },
    };

    inline const Sequences seq_asymm =
    {
        empty,
        empty,
        { "A030228", // n = 2: https://oeis.org/A030228
          { 0, 0, 0, 0, 2, 6, 25, 88, 335, 1215, 4534, 16823, 63159, 237679,
            900341, 3423201, 13073163, 50095285, 192599091, 742576616,
            2870584814, 11122879867, 43191525139, 168046317330, 654998425237,
            2557224396342, 9999083912711, 39153000738695, 153511081627903,
            602621913645490, 2368346964073610, 9317706377210720 } },
        { "A371397", // n = 3: https://oeis.org/A371397
          { 0, 0, 0, 0, 1, 6, 54, 416, 3111, 22898, 168460, 1242985, 9227333,
            68949103, 518618196, 3925228596, 29879207817, 228630283775,
            1757699977107, 13570824097968, 105182547181534, 818093724437992,
            6383353614308209 } },
    };

    inline const Sequences& cubes (const char *name)
    {
        if (! strcmp (name, "free"))  return seq_free;
        if (! strcmp (name, "fixed")) return seq_fixed;
        if (! strcmp (name, "symm"))  return seq_symm;
        if (! strcmp (name, "asymm")) return seq_asymm;
        error ("OEIS cube sequences '%s' not found", name);
    }

    inline const Sequence& cubes (const char *name, int dim)
    {
        const Sequences &ss = cubes (name);
        return dim < (int) ss.size () ? ss[dim] : empty;
    }

    inline Sequence::value_type cubes (const char *name, int dim, int cells)
    {
        return cubes (name, dim)[cells];
    }

    inline Sequence::value_type cubes_fixed (int dim, int cells)
    {
        return cubes ("fixed", dim)[cells];
    }
    inline Sequence::value_type cubes_free (int dim, int cells)
    {
        return cubes ("free", dim)[cells];
    }
    inline Sequence::value_type cubes_symm (int dim, int cells)
    {
        return cubes ("symm", dim)[cells];
    }
    inline Sequence::value_type cubes_asymm (int dim, int cells)
    {
        return cubes ("asymm", dim)[cells];
    }
} // oeis
#endif // OEIS_H

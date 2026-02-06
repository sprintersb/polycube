// -*- c++ -*-
#ifndef PROGRESS_H
#define PROGRESS_H

#include <cstdio>
#include <string>

typedef struct {} progress_with_total_t;
static inline constexpr progress_with_total_t PROGRESS_WITH_TOTAL {};

template<typename T>
struct Progress
{
    static inline constexpr char back[] = "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b"
        "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b";
    T total;
    T margin;
    T last_report = 0;
    bool started = false;
    // 1st % refers to the current progress in T.
    // 2nd % refers to the prgress in percent relative to total.
    std::string fmt;
    Progress (progress_with_total_t, T total, T margin, const char *fmt)
        : total(total), margin(margin), fmt(fmt)
    {
        update ((T) 0);
    }
    Progress (const Progress<T>&) = delete;
    void update (T t)
    {
        if (! started || t - last_report >= margin)
        {
            printf (back);
            printf (fmt.c_str(), t, 100.0 * t / total);
            fflush (stdout);
        }
        started = true;
        last_report = t;
    }
    ~Progress ()
    {
        printf ("%s                                    %s%s",
                back, back, back);
        fflush (stdout);
    }
};
#endif // PROGRESS_H

// -*- c++ -*-
#ifndef PROGRESS_H
#define PROGRESS_H

#include <string>
#include <cstdio>
#include <cassert>

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
    bool done_p = false;
    // 1st % refers to the current progress in T.
    // 2nd % refers to the prgress in percent relative to total.
    std::string fmt;
    Progress (progress_with_total_t, T total, T margin, const char *fmt)
        : total(total), margin(margin), fmt(fmt)
    {
        restart ();
    }
    Progress (const Progress<T>&) = delete;
    void restart ()
    {
        started = done_p = 0;
        update ((T) 0);
    }
    void update (T t)
    {
        assert (! done_p && "done't update a done Progress");
        if (! started || t - last_report >= margin)
        {
            printf (back);
            printf (fmt.c_str(), t, 100.0 * t / total);
            fflush (stdout);
        }
        started = true;
        last_report = t;
    }
    void done ()
    {
        if (! done_p)
        {
            printf ("%s                                    %s%s",
                    back, back, back);
            fflush (stdout);
        }
        done_p = true;
    }
    ~Progress ()
    {
        done ();
    }
};
#endif // PROGRESS_H

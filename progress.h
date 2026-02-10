// -*- c++ -*-
#ifndef PROGRESS_H
#define PROGRESS_H

#include <algorithm> // std::max
#include <string>
#include <sstream>
#include <iostream>
#include <functional>
#include <cassert>
#include <cstdio>

template<typename T>
struct Progress
{
    using value_type = T;
    T total;
    T margin;
    T last_report = 0;
    int len = 0;
    bool started = false;
    bool done_p = false;
    using Printer = std::function<void (std::ostringstream&, Progress&, T)>;
    Printer printer;

    Progress (T total, T margin, Printer printer)
        : total(total), margin(margin), printer(printer)
    {
        reset ();
    }
    Progress (const Progress<T>&) = delete;
    void reset ()
    {
        started = done_p = 0;
        update ((T) 0);
    }
    void update (T t)
    {
        assert (! done_p && "done't update a done Progress");
        if (! started || t - last_report >= margin)
        {
            std::ostringstream stream;
            printer (stream, *this, t);
            erase ();
            len = (int) stream.str().length();
            std::cout << stream.str();
            std::cout.flush ();
        }
        started = true;
        last_report = t;
    }
    static void print_bar (std::ostringstream &ost, int bar_len, double frac)
    {
        char proz[100];
        const int len_text = std::sprintf (proz, " %.2f%%", 100.0 * frac);
        const int len_bars = bar_len - len_text - (int) ost.str().length();
        const int n_bars = std::max<int> (1, frac * len_bars);
        Progress::write (ost, '#', n_bars);
        Progress::write (ost, ' ', len_bars - n_bars);
        ost << proz;
    }
    void done ()
    {
        if (! done_p)
        {
            erase ();
            std::cout.flush ();
        }
        done_p = true;
    }
    static void write (std::ostream &ost, char c, int n)
    {
        for (int i = 0; i < n; ++i)
            ost << c;
    }
    void erase ()
    {
        Progress::write (std::cout, '\b', len);
        Progress::write (std::cout, ' ', len);
        Progress::write (std::cout, '\b', len);
        len = 0;
    }
    ~Progress ()
    {
        done ();
    }
};
#endif // PROGRESS_H

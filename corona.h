// -*- c++ -*-
#ifndef CORONA_H
#define CORONA_H

#include <iostream>
#include <unordered_set>
// Own
#include "dim.h"

struct Corona
{
private:
    using Cells = std::unordered_set<Dim, Dim::Hash>;
    Cells cells;
public:
    int size () const
    {
        return cells.size ();
    }
    void add (Dim d)
    {
        cells.emplace (d);
    }
    void erase (Dim d)
    {
        cells.erase (d);
    }
    bool contains (Dim d) const
    {
        return cells.find (d) != cells.end ();
    }
    struct const_iterator
    {
        Cells::const_iterator it;
        void operator ++ () { ++it; }
        bool operator == (const const_iterator &r) const { return it == r.it; }
        bool operator != (const const_iterator &r) const { return it != r.it; }
        Dim operator * () const { return *it; }
    };
    const_iterator begin () const { return const_iterator { cells.cbegin() }; }
    const_iterator end   () const { return const_iterator { cells.cend() }; }
    const_iterator cbegin () const { return begin (); }
    const_iterator cend   () const { return end (); }
    friend std::ostream& operator << (std::ostream&, const Corona&);
};

inline std::ostream& operator << (std::ostream &ost, const Corona &c)
{
    ost << "{#" << c.size ();
    for (auto c : c.cells)
        ost << " " << c;
    return ost << " }";
}
#endif // CORONA_H

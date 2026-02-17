// -*- c++ -*-
#ifndef ITERATOR_WRAP_H
#define ITERATOR_WRAP_H

#include <cstddef>

template<typename C, typename M, typename S>
struct IteratorWrap
{
    using container_type = C;
    using value_type = typename container_type::value_type;
    using star_type = S;
    using member_type = M;

    IteratorWrap (member_type m) : m_(m) {};

    star_type operator * () const { return *m_; }
    IteratorWrap& operator ++ () { ++m_; return *this; };
    IteratorWrap& operator -- () { --m_; return *this; };
    bool operator == (const IteratorWrap &i) const { return m_ == i.m_; }
    bool operator != (const IteratorWrap &i) const { return m_ != i.m_; }
    bool operator <  (const IteratorWrap &i) const { return m_ < i.m_; }
    ptrdiff_t operator - (const IteratorWrap &i) const { return m_ - i.m_; }
    IteratorWrap operator + (int i) const { return IteratorWrap (m_ + i); }
    IteratorWrap operator - (int i) const { return IteratorWrap (m_ - i); }
private:
    member_type m_;
    IteratorWrap () = delete;
};

#endif // ITERATOR_WRAP_H

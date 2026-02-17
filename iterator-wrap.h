// -*- c++ -*-
#ifndef ITERATOR_WRAP_H
#define ITERATOR_WRAP_H

template<typename C, typename M, typename S>
struct IteratorWrap
{
    using container_type = C;
    using value_type = typename container_type::value_type;
    using star_type = S;
    using member_type = M;
    friend container_type;
    friend value_type;
    void operator ++ () { ++m_; };
    bool operator == (const IteratorWrap &i) const { return m_ == i.m_; }
    bool operator != (const IteratorWrap &i) const { return m_ != i.m_; }
    star_type operator * () const { return *m_; }
    IteratorWrap (member_type m) : m_(m) {};
private:
    member_type m_;
    IteratorWrap () = delete;
};

#endif // ITERATOR_WRAP_H

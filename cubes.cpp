#include <unordered_set>
// C
#include <cstring> // memset
// Own
#include "cubes.h"
#include "debug.h"
#include "util.h"

/////////////////////////////////////////////////////////////////////////////
// Canonicalization

// Mn is used to encode a traversal of the hyperoctahedral group Oh of
// dimension n.  Each letter stands for an elemment of Oh that yields an
// element of Oh not yet encountered.  There is no known formula for such a
// traversal that works in any dimension, hence use the tight encoding below.
// Plus, this is in the hot path, and we want to avoid dups at any cost.
#define M2  "111"
#define M3  M2 "2" M2 "3" M2 "3" M2 "4" M2 "2" M2
#define M4  M3 "5" M3 "5" M3 "6" M3 "5" M3 "6" M3 "5" M3 "5" M3
#define M5  M4 "8" M4 "9" M4 "8" M4 "9" M4 "8" M4 "8" M4 "7" M4 "7" M4 "7" M4

template<typename T>
inline T Cubes::congruents (int dim) const
{
    Assert (dim >= 1 && dim <= DIM, "todo: canonicalization in DIM %d", DIM);
    Cubes c (*this);
    CongruentsAspect<T> aspect;
    aspect.insert (c);
    const char *S = "?";
    switch (dim)
    {
        case 1: S = "0"; break;
        case 2: S = M2 "0" M2; break;
        case 3: S = M3 "0" M3; break;
        case 4: S = M4 "0" M4; break;
        case 5: S = M5 "0" M5; break;
        default:
            unreachable ("todo: canonicalize in DIM %d", DIM);
    }
    for (auto s = S; *s && !aspect.ready(); ++s)
        switch (*s - '0')
        {
            default: unreachable ("todo: implement char '%c' (0x%x)", *s, *s);
            case 0: aspect.insert (c = c.mirrored (0)); break;
            case 1: aspect.insert (c.rotate (0, 1)); break;
#if DIM >= 3
            case 2: aspect.insert (c.rotate (0, 2)); break;
            case 3: aspect.insert (c.rotate (1, 2)); break;
            case 4: aspect.insert (c.rotate (2, 1)); break;
#endif
#if DIM >= 4
            case 5: aspect.insert (c.rotate (3, 0)); break;
            case 6: aspect.insert (c.rotate (3, 1)); break;
#endif
#if DIM >= 5
            case 7: aspect.insert (c.rotate (4, 0)); break;
            case 8: aspect.insert (c.rotate (4, 1)); break;
            case 9: aspect.insert (c.rotate (3, 4)); break;
#endif
        } // switch
    return aspect.value ();
}

// Keeping track of how many cublis have a specific distance to a vertex
// of the bounding box.  The lengths of the diagonals of our bounding boxes
// are all below MAX_DIAGONAL_LENGTH (in Manhattan metric).
// This struct is in the hot path, so we use a managed std::array instead
// of a more dynamic container like std::vector or std::map.
using DistBase = std::array<int16_t, 1 + MAX_DIAGONAL_LENGTH>;
struct Dist : DistBase
{
    int id = -1;
    int value = -1;
    Dist ()
    {
        std::memset (DistBase::data(), 0,
                     DistBase::size() * sizeof (DistBase::value_type));
    }
    int size () const
    {
        return size_;
    }
    value_type& at (int i)
    {
        if (i >= size_)
            resize (1 + i);
        return DistBase::operator [] (i);
    }
    value_type at (int i) const
    {
        return i >= size_ ? 0 : DistBase::operator [] (i);
    }
    value_type& operator [] (int i) { return at (i); }
    value_type  operator [] (int i) const { return at (i); }
    int cmp (const Dist &r) const
    {
        if (size () != r.size ())
            return size () - r.size ();
        for (int i = 0; i < size (); ++i)
            if (at (i) != r.at (i))
                return at (i) - r.at (i);
        return 0;
    }
    bool neighbors_uniquely_sortable (const VertexValues &vvs, int dim) const
    {
        unsigned bits = 0;
        for (int j = 0; j < dim; ++j)
            if (const unsigned mask = 1u << vvs[id ^ (1 << j)]; bits & mask)
                return false;
            else
                bits |= mask;
        return true;
    }
private:
    int size_ = 0;
    void resize (int i)
    {
        Assert (i < (int) DistBase::size (), "%d exceeds MAX_DIAGONAL_LENGTH"
                " %d", (int) DistBase::size (), MAX_DIAGONAL_LENGTH);
        size_ = i;
    }
    iterator begin () = delete;
    iterator end   () = delete;
    const_iterator begin () const = delete;
    const_iterator end   () const = delete;
    const_iterator cbegin () const = delete;
    const_iterator cend   () const = delete;
};

std::ostream& operator << (std::ostream &ost, const Dist &c)
{
    ost << "Dist." << c.id << ":  \t";
    for (int i = 0; i < c.size (); ++i)
        ost << " " << + c.at (i);
    return ost << "\n";
}

// We are only interested in a canonical congruent.
template<>
struct Cubes::CongruentsAspect<Cubes>
{
    Cubes min_cubes;
    Cubes& value () { return min_cubes; }
    const Cubes& value () const { return min_cubes; }
    void insert (const Cubes &c)
    {
        if (min_cubes.size () == 0 || c < min_cubes)
            min_cubes = c;
    }
    bool ready () const { return false; }
};


// We are only interested in the number of congruents.
template<>
struct Cubes::CongruentsAspect<int>
{
    static constexpr int hoho = gjl::hyperoctahedral_order (DIM);
    int size_ = 0;
    int value () const
    {
        // HOh is a group of order hoho, hence when we see more than
        // hoho / 2 elements then we know that it generates all of HOh.
        return size_ > hoho / 2 ? hoho : size_;
    }
    bool ready () const
    {
        return size_ > hoho / 2;
    }
    void insert (const Cubes &c)
    {
        if (! ready ())
            do_insert (c);
    }
private:
    void set_size (int sz)
    {
        Assert (sz <= 1 + hoho / 2, "we need at most 1 + hoho/2 elements");
        size_ = sz;
    }
    // As it seems, even in dimension 5 with hoho = 3840, a managed
    // std::array beats std::unordered_set (and std::set).  Reasons are
    // that std::array lives on the stack and doesn't require expensive
    // heap memory, and it neither requires hashes nor Cubes < Cubes.
#if DIM <= 5
    std::array<Cubes, 1 + hoho / 2> cubs;
    void do_insert (const Cubes &c)
    {
        for (int i = 0; i < size_; ++i)
            if (cubs[i] == c)
                return;
        set_size (1 + size_);
        cubs[size_ - 1] = c;
    }
#else
#warning profile this before using
    struct Hash
    {
        Hasher32::type operator () (const Cubes &c) const
        {
            Hasher32::type h = 0;
            for (Dim d : c)
                h = Hasher32::add (h, d.ival());
            return h;
        }
    };
    std::unordered_set<Cubes, Hash> cubs;
    void do_insert (const Cubes &c)
    {
        cubs.insert (c);
        set_size ((int) cubs.size ());
    }
#endif
};

// Return the dimension number when there is respective simple mirror symmetry.
inline auto Cubes::find_symmetry (const DistPointers &pc, const Box &bbox,
                                  int dim) const -> Symmetry
{
    const Symmetry none; // Just returning {} may raise a GCC warning.
    const int max_value = pc[(1 << dim) - 1]->value;
    if (max_value != (1 << (dim - 1)) - 1)
        return none;
    // # valuations are exacly half of HOh's symmetries, which is a
    // hallmark for mirror symmetry.
    const int symmetry = pc[0]->id ^ pc[1]->id;
    // For now, only consider simple reflections, i.e. symmetry has popcount 1.
    if (symmetry & (symmetry - 1))
        return none;
    for (int i = 2; i < 1 << dim; i += 2)
        if ((pc[i]->id ^ pc[i + 1]->id) != symmetry)
            return none;
    // All looks good until now, but we need a final proof.
    const int d = gjl::count_trailing_zeros (symmetry);
    return matches_flipped (d, bbox) ? d : none;
}

int Cubes::multiplicity () const
{
    Symmetry symmetry;
    if (canonicalizable_vertices (DIM, symmetry))
        return gjl::hyperoctahedral_order (DIM) >> !! symmetry;
    if (DIM >= 3)
    {
        Cubes c (*this);
        Box bbox = c.bounding_box ();
        const int dim = std::max (1, c.squeeze (bbox));
        if (dim == DIM - 1 && c.canonicalizable_vertices (dim, symmetry))
            return DIM * (gjl::hyperoctahedral_order (dim) >> !! symmetry);
    }
    return congruents<int> (DIM);
}

inline bool Cubes::canonicalizable_vertices (int dim, Symmetry &symmetry) const
{
    VertexValues vvs;
    return !! canonical_vertex (vvs, bounding_box (), dim, symmetry);
}

void Cubes::canonicalize ()
{
    Box bbox = bounding_box ();
    Assert (bbox.lo == Dim::all0, "expecting 0-aligned Cubes");
    Assert (bbox.hi.dist (bbox.lo) < (int) std::tuple_size<DistBase>::value,
            "Cubes diameter %d exceeds DistBase capacity of %zd",
            bbox.hi.dist (bbox.lo), std::tuple_size<DistBase>::value);
    const int dim = DIM >= 3
        ? std::max (1, squeeze (bbox))
        : DIM;
    const bool success = maybe_canonicalize_vertices (bbox, dim);
    if (Cubes::take_stat)
    {
        stat[success] += 1;
        stat[2 + success] += success
            ? VERTEX_CANONICALIZATION_COST + (dim != DIM)
            : gjl::hyperoctahedral_order (dim) + (dim != DIM);
    }
    if (! success)
        cells = std::move (congruents<Cubes> (dim).cells);
}

// Return optional id of a canonical vertex in [0, 2^DIM).
// Set optional symmetry to a dimension of mirror symmetry.
inline auto Cubes::canonical_vertex (VertexValues &vvs, const Box &bbox,
                                     int dim, Symmetry &symmetry) const
    -> Vertex
{
    symmetry.reset ();
    // For each cubli, record its Man distances to any of the bounding vertices.
    std::array<Dist, 1 << DIM> vdist;
    for (int id = 0; id < 1 << dim; ++id)
    {
        // The id of a vertex is the canonical integer in [0, 2^DIM).
        Dist &vd = vdist[id];
        vd.id = id;
        const Dim diag = Dim::diag (id);
        const Dim ecke = bbox.vertex (id);
        for (Dim d : *this)
            // Get the distance to bounding vertex id's cube diagonal.
            // This is better than just adding 1 for a cubli of its distance.
            vd[d.dist (ecke)] += 1 + 13 * d.dist (ecke, diag);
    }

    // Sort vdist[].
    DistPointers pc;
    for (int id = 0; id < 1 << DIM /* sic! to avoid warn */; ++id)
        pc[id] = & vdist[id];
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Warray-bounds"
    std::sort (pc.begin (), pc.begin () + (1 << dim),
               [] (const Dist *pa, const Dist *pb) -> bool
               {
                   return pa->cmp (*pb) < 0;
               });
#pragma GCC diagnostic pop
    // Add valuations in increasing order, assigning same valuation
    // if it matches the predecessor's.
    pc[0]->value = 0;
    for (int b = 1; b < 1 << dim; ++b)
        pc[b]->value = pc[b - 1]->value + (pc[b]->cmp (*pc[b - 1]) > 0);
    const int max_value = pc[(1 << dim) - 1]->value;
    // We are after the vertex -> valuation.  Transfer them from auto to output.
    for (int id = 0; id < 1 << dim; ++id)
        vvs[id] = vdist[id].value;
    // Shortcut: When all vertex dists are different, we always succeed.
    if (max_value == (1 << dim) - 1)
        // Canonical: Use the vertex with the smallest valuation.
        return pc[0]->id;

    symmetry = find_symmetry (pc, bbox, dim);
    if (symmetry)
        return pc[0]->id;

    // Now most likely we see a Cubes with some symmetry, but in up
    // to 10% of cases we succeed by looking a bit more closely.
    for (int b = 0; b < 1 << dim; ++b)
    {
        bool good = true;
        good &= b == 0              || pc[b]->value > pc[b - 1]->value;
        good &= b == (1 << dim) - 1 || pc[b]->value < pc[b + 1]->value;
        if (good && pc[b]->neighbors_uniquely_sortable (vvs, dim))
            // Canonical: The vertex with the smallest unique valuation.
            return pc[b]->id;
    }
    return {};
}

inline bool Cubes::maybe_canonicalize_vertices (const Box &bbox, int dim)
{
    VertexValues vvs;
    Symmetry symmetry;
    const Vertex vertex = canonical_vertex (vvs, bbox, dim, symmetry);
    if (! vertex)
        return false;

    bool garbled = false;
    const int id = vertex.value ();
    Assert (vvs.size () == 1 << DIM, "a %dd cube has %d vertices, not %zd",
            DIM, 1 << DIM, vvs.size ());
    Cubes p2 (*this);
    // Reflect such that vertex id becomes vertex 0.
    if (id != 0)
    {
        p2.flip (id, bbox);
        garbled = true;
        // Adjust vvs accordingly.
        for (int b = 0; b < 1 << dim; ++b)
            if (b > (b ^ id))
                std::swap (vvs[b], vvs[b ^ id]);
    }
    // Now sort the dim vertices adjacent to id.
    uint64_t todo = 0;
    for (int b = 0; b < dim; ++b)
        todo |= (uint64_t) 1 << vvs[1 << b];
    // Vertex 0 has dim neighbors, canonicalize them such that higher
    // dimension neighbors get higher vertex values.
    for (int d = 0; d < dim - 1; ++d, todo &= todo - 1)
    {
        Assert (todo, "bad 'todo' popcount");
        const int lsb = gjl::count_trailing_zeros (todo);
        // Search the vertex with the value as specified by todo in
        // the vertex -> value list. We only have to look for vertex
        // indices that are a power of 2 (and hence are a neighbors of 0).
        for (int b = 0; b < dim; ++b)
            if (vvs[1 << b] == lsb
                && b != d)
            {
                // Make vertex b the next neighbor of vertex 0.
                p2.swap (d, b);
                garbled = true;
                // Adjust vvs accordingly.
                for (int i = 0; i < dim; ++i)
                {
                    const int u = 1 << i;
                    const int v = gjl::bitswap (u, d, b);
                    if (u > v)
                        std::swap (vvs[u], vvs[v]);
                }
            }
    }
    // Now p2 has canonical cublis, but swapping and flipping
    // clobbered cells' order.  Re-construct a proper one.
    if (size () <= SORT_THRESHOLD)
    {
        if (garbled)
        {
            cells.clear ();
            for (Dim d : p2)
                add (d);
        }
        else
            cells = std::move (p2.cells);
    }
    else
    {
        cells = std::move (p2.cells);
        if (garbled)
            sort ();
    }
    return true;
}



std::string Cubes::ascii (char c) const
{
#if DIM == 2
    auto bbox = bounding_box ();
    std::string str;
    for (Dim::value_t y = bbox.hi[1]; y >= bbox.lo[1]; --y)
    {
        for (Dim::value_t x = bbox.lo[0]; x <= bbox.hi[0]; ++x)
            str += contains (Dim { x, y }) ? c : ' ';
        str += "\n";
    }
    return str;
#else
    (void) c;
    return "Cubes.ascii(DIM != 2)";
#endif
}

std::ostream& operator << (std::ostream &ost, const Cubes &c)
{
    ost << "{#" << c.size ();
    for (Dim d : c)
        ost << " " << d;
    return ost << " }";
}

#include <unordered_set>
#include <cinttypes>
// Own
#include "cubes.h"
#include "debug.h"
#include "util.h"

#ifdef DEBUG
#define S(x) stat.record_canon (x)
#else
#define S(x) (void) 0
#endif

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
#define M6  M5 "A" M5 "A" M5 "A" M5 "B" M5 "B" M5 "C" M5 "B" M5 "B" \
            M5 "A" M5 "A" M5 "A" M5

template<typename T>
inline T Cubes::congruents (int dim, bool same_parity) const
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
        case 6: S = M6 "0" M6; break;
        default:
            unreachable ("todo: canonicalize in DIM %d", DIM);
    }
    for (auto s = S; *s && !aspect.ready(); ++s)
        switch (*s)
        {
            default: unreachable ("todo: implement char '%c' (0x%x)", *s, *s);
            case '0':
                if (same_parity)
                    return aspect.value ();
                aspect.insert (c = c.mirrored (0));
                break;
            case '1': aspect.insert (c.rotate (0, 1)); break;
#if DIM >= 3
            case '2': aspect.insert (c.rotate (0, 2)); break;
            case '3': aspect.insert (c.rotate (1, 2)); break;
            case '4': aspect.insert (c.rotate (2, 1)); break;
#endif
#if DIM >= 4
            case '5': aspect.insert (c.rotate (3, 0)); break;
            case '6': aspect.insert (c.rotate (3, 1)); break;
#endif
#if DIM >= 5
            case '7': aspect.insert (c.rotate (4, 0)); break;
            case '8': aspect.insert (c.rotate (4, 1)); break;
            case '9': aspect.insert (c.rotate (3, 4)); break;
#endif
#if DIM >= 6
            case 'A': aspect.insert (c.rotate (5, 0)); break;
            case 'B': aspect.insert (c.rotate (5, 2)); break;
            case 'C': aspect.insert (c.rotate (3, 5)); break;
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
        DistBase::fill (0);
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
        if (DIM >= 7)
            unreachable ("todo: DIM >= 7");
        uint64_t bits = 0;
        for (int j = 0; j < dim; ++j)
            if (uint64_t mask = (uint64_t) 1 << vvs[id ^ (1 << j)]; bits & mask)
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
    // As it seems, even in dimension 4 with hoho = 384, a managed
    // std::array beats std::unordered_set (and std::set).  Reasons are
    // that std::array lives on the stack and doesn't require expensive
    // heap memory, and it neither requires hashes nor Cubes < Cubes.
#if DIM <= 4
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

// Return 2 (mirror symmetry) or 0 (unknown symmetry).
inline auto Cubes::find_symmetry (const DistPointers &pc, const Box &bbox,
                                  int dim) const -> Symmetry
{
    const int max_value = pc[(1 << dim) - 1]->value;
    if (max_value != (1 << (dim - 1)) - 1)
        return 0;
    // # valuations are exacly half of HOh's symmetries, which is a
    // hallmark for mirror symmetry.
    const int symmetry = pc[0]->id ^ pc[1]->id;
    // For now, only consider simple reflections, i.e. symmetry has popcount 1.
    if (symmetry & (symmetry - 1))
        return 0;
    for (int i = 2; i < 1 << dim; i += 2)
        if ((pc[i]->id ^ pc[i + 1]->id) != symmetry)
            return 0;
    // All looks good until now, but we need a final proof.
    const int d = gjl::count_trailing_zeros (symmetry);
    return matches_flipped (d, bbox) ? 2 : 0;
}

// Return 2 (diagonal symmetry) or 0 (unknown symmetry).
inline bool Cubes::is_diagonal_symmetric (int id1, int id2,
                                          const Context &ctx) const
{
    const int ix = id1 ^ id2;
    // Avoiding popcount.
    const int ix1 = ix & (ix - 1); // ix without LSB.
    const int ix2 = ix1 & (ix1 - 1); // ix without 2 LSBs.
    // id1 and id2 should differ in exactly 2 bits for this kind of symmetry:
    // Fail if popcount (ix) != 2.
    if (ix1 == 0 || ix2 != 0)
        return false;
    // id1 and id2 differ in two bits, and we have two symmetry cases:
    // *0*1* <-> *1*0* = reflection at main diagonal y = x.
    // *0*0* <-> *1*1* = reflection at minor diagonal y = hi - x.
    const int a = gjl::count_trailing_zeros (ix);
    const int b = gjl::count_trailing_zeros (ix1);
    const int hi = ctx.bbox.hi[a];
    if (hi == ctx.bbox.hi[b])
        return id2 == (id1 | ix) // Notice that id2 > id1.
            // Reflect at minor diagonal.
            ? matches_flipped_swapped (a, b, hi)
            // Reflect at main diagonal: Swap coordinates a <-> b.
            : matches_swapped (a, b);
    return false;
}

inline auto Cubes::find_diag_symmetry (const DistPointers &pc,
                                       const Context &ctx) const -> Symmetry
{
    const int dim = ctx.dim;
    const int max_value = pc[(1 << dim) - 1]->value;
    // Only 75% of possible colors is a hallmark of diagonal symmetry.
    // These are the diagonal symmetries of the square together with
    // higher dimensional versions thereof.
    if (dim >= 2 && max_value == (1 << dim) - (1 << dim) / 4 - 1)
        // Find 2 vertices of same color.
        for (int i = 1; i < 1 << dim; ++i)
            if (pc[i]->value == pc[i - 1]->value)
            {
                // Any vertex pair with like colors determines the symmetry.
                const int id1 = std::min (pc[i - 1]->id, pc[i]->id);
                const int id2 = std::max (pc[i - 1]->id, pc[i]->id);
                return is_diagonal_symmetric (id1, id2, ctx) ? 2 : 0;
            }
    return 0;
}

// Size of the PolyCube's orbit in HOh.  Symmetry is known to be:
// 2) Mirror symmetric, i.e. multiplicity <= #HOh / 2,
// 1) Asymmetric, i.e. multiplicity = #HOh,
// 0) Not mirror symmetric; more is not known.
int Cubes::multiplicity (Symmetry symmetry) const
{
    if (symmetry == 1)
        return gjl::hyperoctahedral_order (DIM);
    bool symmetric = symmetry == 2;
    Context ctx;
    ctx.dim = DIM;
    ctx.bbox = bounding_box ();
    if (canonical_vertex (ctx))
        return gjl::hyperoctahedral_order (DIM) / ctx.symmetry;
    if (DIM >= 3)
    {
        Cubes c (*this);
        ctx.bbox = c.bounding_box ();
        ctx.dim = std::max (1, c.squeeze (ctx.bbox));
        if (ctx.dim == DIM - 1 && c.canonical_vertex (ctx))
            return DIM * (gjl::hyperoctahedral_order (ctx.dim) / ctx.symmetry);
        if (ctx.dim < DIM)
            symmetric = true;
    }
    // Known mirror symmetry improves speed.
    return congruents<int> (DIM, symmetric);
}

void Cubes::canonicalize ()
{
    Context ctx;
    ctx.bbox = bounding_box ();
    Assert (ctx.bbox.lo == Dim::all0, "expecting 0-aligned Cubes");
    [[maybe_unused]] const int diag_length = ctx.bbox.hi.dist (ctx.bbox.lo);
    Assert (diag_length <= MAX_DIAGONAL_LENGTH, "Cubes diameter %d exceeds "
            "DistBase capacity of %d", diag_length, MAX_DIAGONAL_LENGTH);
    ctx.dim = DIM >= 3
        ? std::max (1, squeeze (ctx.bbox))
        : DIM;
    const bool success = maybe_canonicalize_vertices (ctx);
    if (! success)
        cells = std::move (congruents<Cubes> (ctx.dim, 0).cells);

    if (Cubes::take_stat)
        stat.record_success (success, ctx.dim);
}

// Find a canonical representation, and determine whether it is congruent
// to some flipped version of itself.
// Symmetry = 2: Has mirror symmetry, i.e. multiplicity <= #HOh / 2.
// Symmetry = 1: Has no symmetry, i.e. multiplicity = #HOh.
// Symmetry = 0: Doesn't have mirror symmetry.
void Cubes::canonicalize (Symmetry &symmetry)
{
    Context ctx;
    ctx.bbox = bounding_box ();
    ctx.dim = DIM >= 3
        ? std::max (1, squeeze (ctx.bbox))
        : DIM;
    const bool success = maybe_canonicalize_vertices (ctx);
    if (success)
    {
        symmetry = ctx.dim < DIM || ctx.symmetry == 2 ? 2 : 1;
    }
    // When vertex canonicalization fails then go brute force.
    else if (ctx.dim < DIM)
    {
        cells = std::move (congruents<Cubes> (ctx.dim, 0).cells);
        symmetry = 2;
    }
    else
    {
        cells = std::move (congruents<Cubes> (DIM, 1).cells);
        Cubes &&c2 = mirrored (0).congruents<Cubes> (DIM, 1);
        const int j = cmp (c2);
        symmetry = j == 0 ? 2 : 0;
        if (j > 0)
            cells = std::move (c2.cells);
    }

    if (Cubes::take_stat)
        stat.record_success (success, ctx.dim);
}

// Return the smallest index in pc[] (if any) that points to a unique color.
inline int Cubes::uniquely_valued (const DistPointers &pc, int dim, bool flt,
                                   const Context &ctx)
{
    for (int b = 0; b < 1 << dim; ++b)
        if (b == 0 || pc[b]->value > pc[b - 1]->value)
            if (b == (1 << dim) - 1 || pc[b]->value < pc[b + 1]->value)
                if (! flt || pc[b]->neighbors_uniquely_sortable (ctx.vvs, dim))
                    return b;
    return -1;
}

// On success:
// 1) Set ctx.vertex to the id of a canonical vertex in [0, 2^DIM),
// 2) Set ctx.vvs[id] to the values / colors assigned to the 2^DIM vertices,
// 3) Set ctx.symmetry to 1 or 2: We have multiplicity #HOh / ctx.symmetry.
// On failure, ctx.vertex is empty.  In any case ctx.bbox is unchanged.
inline bool Cubes::canonical_vertex (Context &ctx) const
{
    const int dim = ctx.dim;
    Assert (ctx.bbox == bounding_box (), "bad box");

    // For each cubli, record its distances to any of the bounding vertices.
    VertexTraits<Dist> vdist;
    for (int id = 0; id < 1 << dim; ++id)
    {
        // The id of a vertex is the canonical integer in [0, 2^DIM).
        Dist &vd = vdist[id];
        vd.id = id;
        const Dim diag = Dim::diag (id);
        const Dim ecke = ctx.bbox.vertex (id);
        for (Dim d : *this)
            // Get the distance to bounding vertex id's cube diagonal.
            // This is better than just adding 1 for a cubli of its distance.
            //vd[d.dist (ecke)] += 1 + 13 * d.dist (ecke, diag);
            vd[d.dist (ecke)] += 1 + 13 * d.pseudo_dist (ecke, diag);
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
        ctx.vvs[id] = vdist[id].value;
    // Shortcut: When all vertex dists are different, we always succeed.
    if (max_value == (1 << dim) - 1)
        // Canonical: Use the vertex with the smallest valuation.
        return S(1), ctx.vertex = pc[0]->id, ctx.symmetry = 1, true;

    if (find_symmetry (pc, ctx.bbox, dim))
        return S(2), ctx.vertex = pc[0]->id, ctx.symmetry = 2, true;

    if (find_diag_symmetry (pc, ctx))
    {
        // Canonical: Use the smallest vertex with a unique color.
        // Neighbors are uniquely sortable (modulo diag symmetry).
        const int b = uniquely_valued (pc, dim, false, ctx);
        Assert (b >= 0, "id must exist");
        return S(3), ctx.vertex = pc[b]->id, ctx.symmetry = 2, true;
    }

    // Now most likely we see a Cubes with some symmetry, but in up
    // to 10% of cases we succeed by looking a bit more closely.
    const int b = uniquely_valued (pc, dim, true, ctx);
    if (b >= 0)
        // Canonical: The vertex with the smallest unique valuation.
        return S(4), ctx.vertex = pc[b]->id, ctx.symmetry = 1, true;

    return S(0), ctx.vertex = Vertex {}, ctx.symmetry = 0, false;
}

inline bool Cubes::maybe_canonicalize_vertices (Context &ctx, bool same_parity)
{
    stat.take_canon = 1;
    if (canonical_vertex (ctx))
        canonicalize_vertices (ctx, same_parity);
    stat.take_canon = 0;
    return !! ctx.vertex;
}

// Carry out vertex canonicalization by their coloring in ctx.vvs[].
inline void Cubes::canonicalize_vertices (Context &ctx, bool same_parity)
{
    Assert (!! ctx.vertex, "empty Vertex");
    Assert (ctx.bbox == bounding_box (), "bad box");
    const int id = *ctx.vertex;
    const int dim = ctx.dim;
    VertexValues &vvs = ctx.vvs;

    bool garbled = false;
    bool parity = 0;
    Assert (vvs.size () == 1 << DIM, "a %dd cube has %d vertices, not %zd",
            DIM, 1 << DIM, vvs.size ());
    Cubes p2 (*this);
    // Reflect such that vertex id becomes vertex 0.
    if (id != 0)
    {
        p2.flip (id, ctx.bbox);
        parity = gjl::popcount (id);
        garbled = true;
        // Adjust vvs accordingly.
        for (int b = 0; b < 1 << dim; ++b)
            if (b > (b ^ id))
                std::swap (vvs[b], vvs[b ^ id]);
    }
    // Now sort the dim vertices adjacent to id.
    uint64_t todo = 0;
    if (DIM >= 7)
        unreachable ("todo: DIM >= 7");
    for (int b = 0; b < dim; ++b)
        todo |= (uint64_t) 1 << vvs[1 << b];
    // Vertex 0 has dim neighbors, canonicalize them such that higher
    // dimension neighbors get higher vertex values.
    // !!! Notice that this works as expected even in the case of main
    // !!! diagonal symmetry where exactly two neighbors have like colors.
    for (int d = 0; d < dim - 1; ++d, todo &= todo - 1)
    {
        Assert (todo, "bad 'todo' popcount");
        const int lsb = gjl::count_trailing_zeros (todo);
        // Search the vertex with the value as specified by todo in
        // the vertex -> value list. We only have to look for vertex
        // indices that are a power of 2 (and hence are neighbors of 0).
        for (int b = 0; b < dim; ++b)
            if (vvs[1 << b] == lsb
                && b != d)
            {
                // Make vertex b the next neighbor of vertex 0.
                p2.swap (d, b);
                parity ^= 1;
                garbled = true;
                // Adjust vvs accordingly.
                for (int i = 0; i < dim; ++i)
                {
                    const int u = 1 << i;
                    const int v = gjl::bitswap (u, d, b);
                    if (u > v)
                        std::swap (vvs[u], vvs[v]);
                }
                break;
            }
    }
    if (same_parity && parity)
    {
        p2.flip (1 << 0);
        // Adjust vvs accordingly for future.
        for (int b = 0; b < 1 << dim; b += 2)
            std::swap (vvs[b], vvs[b ^ 1]);
    }
    // Now p2 has canonical cublis, but swapping and flipping
    // clobbered cublis' order.  Re-construct a proper one.
    if (size () < SORT_THRESHOLD)
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


void Cubes::Stat::record_success (bool success, int dim)
{
    n_succ[success] += 1;
    n_succ[2 + success] += success
        ? VERTEX_CANONICALIZATION_COST + (dim != DIM)
        : gjl::hyperoctahedral_order (dim) + (dim != DIM);
}

inline void Cubes::Stat::record_canon (int id)
{
    n_canon.at (id) += 1;
}

void Cubes::Stat::print () const
{
    enum { STAT_FAIL, STAT_SUCC, STAT_COST_FAIL, STAT_COST_SUCC };
    // Stat about fraction of fast canonicalization.
    const int64_t s0 = n_succ[STAT_FAIL];
    const int64_t s1 = n_succ[STAT_SUCC];
    const int64_t s2 = n_succ[STAT_COST_FAIL];
    const int64_t s3 = n_succ[STAT_COST_SUCC];
    auto tot = s0 + s1;
    if (tot)
    {
        const double f0 = tot ? (double) s0 / tot : 0.0;
        const double f1 = tot ? (double) s1 / tot : 0.0;
        printf ("Stat: %" PRIi64 "(%.2f%%)  %" PRIi64 "(%.2f%%)",
                s0, 100 * f0, s1, 100 * f1);
        const double cf0 = s0 ? (double) s2 / s0 : 0.0;
        const double cf1 = s1 ? (double) s3 / s1 : 0.0;
        printf ("\t Cost factor: %.1f + %.1f = %.1f\n",
                cf0, cf1, f0 * cf0 + f1 * cf1);
    }
    else
        printf ("Cost factor: %.1f\n", 0. + gjl::hyperoctahedral_order (DIM));

    int64_t ss = 0;
    for (auto &a : n_canon)
        ss += a;
    if (ss)
    {
        std::vector<double> dd (n_canon.size());
        for (size_t i = 0; i < n_canon.size (); ++i)
            dd[i] = 100.0 * n_canon[i] / ss;
        printf ("Cano: %.1f%% / %.1f%% + %.1f%% / %.1f%% + %.1f%%\n",
                dd[0], dd[1], dd[4], dd[2], dd[3]);
    }
}

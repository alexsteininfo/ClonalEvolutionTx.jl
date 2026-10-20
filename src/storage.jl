# Saving and loading populations in a small, stable binary format, and Newick export.
#
# A file is a header followed by raw little-endian columns (every array goes through
# `htol` on writing and `ltoh` on reading); it does not depend on Julia's Serialization.
# Layout, version 1 (all integers little-endian):
#   magic             8 bytes  "CEpop\0\0\0"
#   version           UInt32   1
#   mode              UInt8    1 tree, 2 cells
#   t, m              Float64
#   clonal_mutations  Int64
#   next_id           Int64
#   tree mode:  n Int64, then parent, left, right (Int32), id (Int64), birthtime (Float64),
#               clone, neutral, total (Int32), each n long, without removed rows;
#               then na Int64 and alive (Int32, na long)
#   cell mode:  n Int64, then id (Int64), birthtime (Float64), clone, divisions, neutral
#               (Int32), each n long
#   clone table: nk Int64, then parent (Int32), s, α (Float64), kind (UInt8: 1 founder,
#               2 driver, 3 resistance, 4 custom), origin_time (Float64), origin_cell,
#               size (Int64), each nk long

const _MAGIC = UInt8['C', 'E', 'p', 'o', 'p', '\0', '\0', '\0']
const _FORMAT_VERSION = UInt32(1)
const _KINDS = (:founder, :driver, :resistance, :custom)

_kind_code(k::Symbol) = UInt8(something(findfirst(==(k), _KINDS)))

_writecol(io::IO, v::AbstractVector) = write(io, htol.(v))
_writecount(io::IO, n::Integer) = write(io, htol(Int64(n)))

function _write_clones(io::IO, ct::CloneTable)
    _writecount(io, n_clones(ct))
    _writecol(io, ct.parent)
    _writecol(io, ct.s)
    _writecol(io, ct.α)
    _writecol(io, UInt8[_kind_code(k) for k in ct.kind])
    _writecol(io, ct.origin_time)
    _writecol(io, ct.origin_cell)
    _writecol(io, ct.size)
end

function _write_header(io::IO, mode::UInt8, pop::AbstractPopulation)
    write(io, _MAGIC)
    write(io, htol(_FORMAT_VERSION))
    write(io, mode)
    write(io, htol(Float64(pop.t)))
    write(io, htol(Float64(pop.m)))
    write(io, htol(Int64(pop.clonal_mutations)))
    write(io, htol(Int64(pop.next_id)))
end

"""
    save_population(path, pop::AbstractPopulation) -> path

Write `pop` to `path` in the package's binary format and return `path`. The file holds,
in little-endian raw columns, a header (magic `"CEpop\\0\\0\\0"`, `UInt32` format version
`1`, `UInt8` mode, `t`, `m`, `clonal_mutations`, `next_id`), then the cells (tree mode:
the eight lineage-tree columns without removed rows, then the living-cell indices; cell
mode: the five cell columns) and the clone table (`parent`, `s`, `α`, `kind` as `UInt8`
codes 1–4, `origin_time`, `origin_cell`, `size`). The neutral stream and the log are not
stored, and `pop` itself is never modified (a tree is compacted in a copy). Read it back
with [`load_population`](@ref).
"""
function save_population(path::AbstractString, pop::TreePopulation)
    tree, alive = pop.tree, pop.alive
    if tree.ndead > 0
        tree = copy(tree)
        remap = _compact!(tree, Int32(1))
        alive = Int32[remap[i] for i in alive]
    end
    open(path, "w") do io
        _write_header(io, 0x01, pop)
        _writecount(io, length(tree))
        for v in (tree.parent, tree.left, tree.right, tree.id, tree.birthtime, tree.clone,
                  tree.neutral, tree.total)
            _writecol(io, v)
        end
        _writecount(io, length(alive))
        _writecol(io, alive)
        _write_clones(io, pop.clones)
    end
    return path
end

function save_population(path::AbstractString, pop::CellPopulation)
    cells = pop.cells
    open(path, "w") do io
        _write_header(io, 0x02, pop)
        _writecount(io, length(cells))
        _writecol(io, Int64[c.id for c in cells])
        _writecol(io, Float64[c.birthtime for c in cells])
        _writecol(io, Int32[c.clone for c in cells])
        _writecol(io, Int32[c.divisions for c in cells])
        _writecol(io, Int32[c.neutral for c in cells])
        _write_clones(io, pop.clones)
    end
    return path
end

"""
    save_tree(path, pop::TreePopulation, node::CellNode) -> path

Write the subtree under `node` (a cell of `pop.tree`) to `path` as a single-root
tree-mode population file (the format of [`save_population`](@ref)). `node` becomes the
root and carries every neutral mutation on its path from the founder: in the file its
`neutral` count is set to its `total`, so the mutations inherited from above count at
`k = n` in the site frequency spectrum and are part of `clonal_mutations` and of every
cell's `mutations_per_cell`, which equal the source cells' burdens. `divisions_per_cell`
of the loaded file counts the divisions from the subtree root. The living cells of the
file are the living cells under `node`, in the order of `pop.alive`; the clone table is
copied with sizes recounted over those cells; `t`, `m`, `clonal_mutations` and `next_id`
are those of `pop`. A `node` that does not belong to `pop.tree`, or whose lineage has been
removed, gives an `ArgumentError`.
"""
function save_tree(path::AbstractString, pop::TreePopulation, node::CellNode)
    tree = pop.tree
    node.tree === tree || throw(ArgumentError("the node does not belong to the population's tree"))
    r = _index(node)
    _isdead(tree, r) && throw(ArgumentError("cell $(node.id) was removed from the tree"))
    nodes, _ = _preorder_idx(tree, r)
    kept = sort!(nodes)
    newtree = _copy_rows(tree, kept)
    newtree.neutral[1] = newtree.total[1]     # row 1 is `node`: its descendants come after it
    insub = falses(length(tree))
    insub[kept] .= true
    alive = Int32[searchsortedfirst(kept, i) for i in pop.alive if insub[i]]
    clones = _recounted_clones(pop.clones, Int32[newtree.clone[i] for i in alive])
    sub = TreePopulation(newtree, alive, clones, pop.t, pop.m, pop.clonal_mutations,
                         pop.next_id, Xoshiro(0), BlockSummary[])
    return save_population(path, sub)
end

# ── Loading ───────────────────────────────────────────────────────────────────

function _read_clones(io::IO, path, readcol, readcount)
    nk = readcount()
    parent = readcol(Int32, nk)
    s = readcol(Float64, nk)
    α = readcol(Float64, nk)
    codes = readcol(UInt8, nk)
    origin_time = readcol(Float64, nk)
    origin_cell = readcol(Int64, nk)
    size = readcol(Int64, nk)
    all(c -> 1 <= c <= length(_KINDS), codes) ||
        throw(ArgumentError("$path: unknown clone kind code"))
    ct = CloneTable()
    for k in 1:nk
        0 <= parent[k] < k || throw(ArgumentError("$path: invalid clone parent for clone $k"))
        size[k] >= 0 || throw(ArgumentError("$path: negative size for clone $k"))
        β, δ = _checked_rates(ct.birth, ct.death, s[k], k)
        push!(ct.β, β); push!(ct.δ, δ)
    end
    append!(ct.parent, parent); append!(ct.s, s); append!(ct.α, α)
    append!(ct.kind, Symbol[_KINDS[c] for c in codes])
    append!(ct.origin_time, origin_time); append!(ct.origin_cell, origin_cell)
    append!(ct.size, size)
    ct.alive_pos = zeros(Int32, nk)
    ct.ever_alive = fill(true, nk)       # no loaded size-0 clone may be revived
    return _rebuild_derived!(ct)
end

"""
    load_population(path) -> AbstractPopulation

Read a file written by [`save_population`](@ref) or [`save_tree`](@ref): a
[`TreePopulation`](@ref) or [`CellPopulation`](@ref) according to the file's mode. A file
with the wrong magic, an unknown format version, mode or clone-kind code, inconsistent
indices, or that ends early gives an `ArgumentError` naming the file.

The clone table is rebuilt from the file with its default rates (`β ≡ 1`, `δ ≡ 0`) and every
clone counts as having lived, so a clone with no cells cannot be revived; the next
`simulate!` block rebinds the rates. The neutral stream is not stored: the loaded population
gets `Xoshiro(0)`, which the next block reseeds, and an empty log. Hence a run continued from
a loaded population is not bit-identical to the uninterrupted run (SPEC §8), though it is
statistically equivalent.
"""
function load_population(path::AbstractString)
    open(path, "r") do io
        try
            return _read_population(io, path)
        catch e
            e isa EOFError && throw(ArgumentError("$path is truncated: the file ended early"))
            rethrow()
        end
    end
end

function _read_population(io::IO, path::AbstractString)
    total = filesize(path)
    magic = read(io, length(_MAGIC))
    magic == _MAGIC || throw(ArgumentError("$path is not a ClonalEvolutionTx population file"))
    version = ltoh(read(io, UInt32))
    version == _FORMAT_VERSION || throw(ArgumentError(
        "$path has format version $version; this package reads version $_FORMAT_VERSION"))
    mode = read(io, UInt8)
    (mode == 0x01 || mode == 0x02) || throw(ArgumentError("$path has unknown mode code $mode"))
    t = ltoh(read(io, Float64))
    m = ltoh(read(io, Float64))
    clonal = Int(ltoh(read(io, Int64)))
    next_id = ltoh(read(io, Int64))
    readcount() = begin
        len = ltoh(read(io, Int64))
        len >= 0 || throw(ArgumentError("$path has a negative length"))
        Int(len)
    end
    readcol(T, len) = begin
        len * sizeof(T) <= total - position(io) && return ltoh.(read!(io, Vector{T}(undef, len)))
        throw(ArgumentError("$path is truncated: the file ended early"))
    end
    n = readcount()
    if mode == 0x01
        parent, left, right = readcol(Int32, n), readcol(Int32, n), readcol(Int32, n)
        id, birthtime = readcol(Int64, n), readcol(Float64, n)
        clone, neutral, tot = readcol(Int32, n), readcol(Int32, n), readcol(Int32, n)
        alive = readcol(Int32, readcount())
        ct = _read_clones(io, path, readcol, readcount)
        inrange(v, lo) = all(x -> lo <= x <= n, v)
        (inrange(parent, 0) && inrange(left, 0) && inrange(right, 0) && inrange(alive, 1)) ||
            throw(ArgumentError("$path: a tree index is out of range"))
        all(k -> 1 <= k <= n_clones(ct), clone) ||
            throw(ArgumentError("$path: a cell has an unknown clone"))
        tree = LineageTree(parent, left, right, id, birthtime, clone, neutral, tot, 0)
        return TreePopulation(tree, alive, ct, t, m, clonal, next_id, Xoshiro(0), BlockSummary[])
    else
        id, birthtime = readcol(Int64, n), readcol(Float64, n)
        clone, divisions, neutral = readcol(Int32, n), readcol(Int32, n), readcol(Int32, n)
        ct = _read_clones(io, path, readcol, readcount)
        all(k -> 1 <= k <= n_clones(ct), clone) ||
            throw(ArgumentError("$path: a cell has an unknown clone"))
        cells = Cell[Cell(id[i], birthtime[i], clone[i], divisions[i], neutral[i]) for i in 1:n]
        return CellPopulation(cells, ct, t, m, clonal, next_id, Xoshiro(0), BlockSummary[])
    end
end

# ── Newick ────────────────────────────────────────────────────────────────────

# Follow the single-child chain from `i`; returns the first node that is a leaf or has two
# children, and the number of edges walked.
function _collapse(tree::LineageTree, i::Int32)
    steps = 0
    @inbounds while true
        l, r = tree.left[i], tree.right[i]
        (l == NOPARENT) == (r == NOPARENT) && return i, steps   # leaf or branching
        i = l == NOPARENT ? r : l
        steps += 1
    end
end

# The time of a collapsed node on the Newick time axis: a leaf (a living cell) sits at the
# population's time `t`, a branching node at its division, the birth time of its daughters.
@inline _newick_time(tree::LineageTree, i::Int32, t::Float64) =
    @inbounds _isleaf(tree, i) ? t : tree.birthtime[tree.left[i]]

"""
    write_newick(io::IO, pop::TreePopulation; branch_lengths = :divisions)
    write_newick(path::AbstractString, pop::TreePopulation; branch_lengths = :divisions)

Write the genealogy of the living cells as Newick, one tree per root (founder) and one
line per tree, each ending in `;` and a newline (an empty population writes nothing).
Leaves are labelled by cell id; chains of single-child nodes are collapsed. Each tree is
rooted at the most recent common ancestor of the cells below its founder, and the root's
branch length is its stem, from the founder to that ancestor. `branch_lengths` is
`:divisions` or `:time` (anything else is an `ArgumentError`):

- `:divisions`: a branch length is the number of divisions along the collapsed chain, and
  the stem the divisions from the founder to the ancestor, e.g. `(4:2,8:3):0;` or
  `(8:1,9:1):2;`.
- `:time`: every node sits on the time axis, a branching node at its division (the birth
  time of its daughters) and a leaf at `pop.t`, so a branch length is the time between a
  node and the division above it, and the stem is the ancestor's division time minus the
  founder's birth time. The tree is ultrametric (every root-to-leaf length is `pop.t`
  minus the founder's birth time), and the length from a leaf up to the common ancestor of
  two cells is their [`coalescence_times`](@ref) entry.

A founder with a single living cell gives `(17:0):d;`, where `d` is its stem: the
divisions from the founder to the cell, or `pop.t` minus the founder's birth time. Written
iteratively, so deep trees are safe.
"""
function write_newick(io::IO, pop::TreePopulation; branch_lengths::Symbol = :divisions)
    branch_lengths in (:divisions, :time) || throw(ArgumentError(
        "branch_lengths must be :divisions or :time, got :$branch_lengths"))
    bytime = branch_lengths === :time
    tree, t = pop.tree, pop.t
    _fmt(x) = bytime ? string(x) : string(Int(x))
    # Stack items: (0, node, length) a collapsed node to write; (1, _, length) the closing
    # parenthesis with the node's branch length; (2, _, _) a comma.
    stack = Tuple{Int8,Int32,Float64}[]
    for r in _rootidxs(tree)
        top, steps = _collapse(tree, r)
        stem = bytime ? _newick_time(tree, top, t) - tree.birthtime[r] : Float64(steps)
        if _isleaf(tree, top)
            print(io, "(", tree.id[top], ":0):", _fmt(stem), ";\n")
            continue
        end
        push!(stack, (Int8(0), top, stem))
        while !isempty(stack)
            kind, i, ℓ = pop!(stack)
            if kind == 2
                print(io, ",")
            elseif kind == 1
                print(io, "):", _fmt(ℓ))
            elseif _isleaf(tree, i)
                print(io, tree.id[i], ":", _fmt(ℓ))
            else
                print(io, "(")
                push!(stack, (Int8(1), Int32(0), ℓ))
                kids = Tuple{Int32,Float64}[]
                ti = bytime ? _newick_time(tree, i, t) : 0.0
                for c in (tree.left[i], tree.right[i])
                    j, s = _collapse(tree, c)
                    push!(kids, (j, bytime ? _newick_time(tree, j, t) - ti : Float64(s + 1)))
                end
                push!(stack, (Int8(0), kids[2][1], kids[2][2]))
                push!(stack, (Int8(2), Int32(0), 0.0))
                push!(stack, (Int8(0), kids[1][1], kids[1][2]))
            end
        end
        print(io, ";\n")
    end
    return nothing
end

function write_newick(path::AbstractString, pop::TreePopulation; branch_lengths::Symbol = :divisions)
    open(path, "w") do io
        write_newick(io, pop; branch_lengths)
    end
    return nothing
end

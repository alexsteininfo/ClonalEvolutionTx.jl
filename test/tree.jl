# The struct-of-arrays lineage tree: handles, pruning, compaction, traversals.

const CE = ClonalEvolutionTx
const P = CE._push_node!

# 1 (root) -> 2 (left) -> {4, 5}; 1 -> 3 (right). Leaves in order: 4, 5, 3.
function tree_fixture()
    tree = LineageTree()
    r = P(tree, CE.NOPARENT, 1, 0.0, 1, 5, 5)
    a = P(tree, r, 2, 1.0, 1, 1, 6); tree.left[r] = a
    b = P(tree, r, 3, 1.0, 2, 2, 7); tree.right[r] = b
    c = P(tree, a, 4, 2.0, 1, 0, 6); tree.left[a] = c
    d = P(tree, a, 5, 2.0, 3, 3, 9); tree.right[a] = d
    return tree
end

@testset "push, links and path totals" begin
    tree = LineageTree()
    r = CE._push_node!(tree, CE.NOPARENT, 1, 0.0, 1, 2, 2)
    a = CE._push_node!(tree, r, 2, 1.0, 1, 1, 3); tree.left[r] = a
    b = CE._push_node!(tree, r, 3, 1.0, 2, 3, 5); tree.right[r] = b
    @test CE._leaves_idx(tree, r) == [a, b]
    @test CellNode(tree, 3).clone == 2 && CellNode(tree, 3).total == 5
    @test length(collect(Leaves(CellNode(tree, 1)))) == 2
end

@testset "CellNode properties read the tree's columns" begin
    tree = tree_fixture()
    root = CellNode(tree, 1)
    @test root.id == 1 && root.birthtime == 0.0 && root.clone == 1
    @test root.neutral == 5 && root.total == 5
    @test isnothing(root.parent)
    L = root.left
    @test L.id == 2 && L.parent == root && L.total == 6
    @test root.right.id == 3 && isnothing(root.right.left) && root.right.clone == 2
    @test [l.id for l in Leaves(root)] == [4, 5, 3]
    @test [n.id for n in PreOrderDFS(root)] == [1, 2, 4, 5, 3]
    @test [c.id for c in children(root)] == [2, 3]
    @test nodevalue(L) == CE.TreeCell(2, 1.0, 1, 1, 6)
    @test_throws ArgumentError root.bogus
    @test sprint(print_tree, root) isa String
    @test length(tree) == 5
end

@testset "equality and hashing are by tree and id" begin
    tree = tree_fixture()
    root = CellNode(tree, 1)
    @test root.left.left == root.left.left
    @test hash(root.left) == hash(root.left)
    @test root.left != root.right
    @test root != CellNode(tree_fixture(), 1)
    @test length(Set([root, root.left, root.left, root.right])) == 3
end

@testset "index traversals" begin
    tree = tree_fixture()
    nodes, ppos = CE._preorder_idx(tree, Int32(1))
    @test nodes == [1, 2, 4, 5, 3] && ppos == [0, 1, 2, 2, 1]
    @test CE._leafcounts(tree, nodes, ppos) == [3, 2, 1, 1, 1]
    @test CE._leaves_idx(tree, Int32(2)) == [4, 5]
    @test CE._rootidx(tree, 5) == 1
    @test CE._rootidxs(tree) == [1]
    @test CE._alive_idx(tree) == [3, 4, 5]
end

@testset "prune keeps unary nodes, compact keeps order" begin
    tree = LineageTree()
    r = P(tree, CE.NOPARENT, 1, 0.0, 1, 0, 0)
    a = P(tree, r, 2, 1.0, 1, 0, 0); tree.left[r] = a
    b = P(tree, r, 3, 1.0, 1, 0, 0); tree.right[r] = b
    c = P(tree, a, 4, 2.0, 1, 0, 0); tree.left[a] = c
    d = P(tree, a, 5, 2.0, 1, 0, 0); tree.right[a] = d
    CE._prune!(tree, b); CE._prune!(tree, c)
    @test CE._isdead(tree, b) && CE._isdead(tree, c)
    @test !CE._isdead(tree, r) && !CE._isdead(tree, a)   # unary nodes stay
    @test CE._leaves_idx(tree, r) == [d]
    @test tree.ndead == 2
    remap = CE._compact!(tree, Int32(1))
    @test tree.id == [1, 2, 5] && remap[d] == 3 && tree.ndead == 0
    @test tree.parent == [0, 1, 2] && tree.left == [2, 0, 0] && tree.right == [0, 3, 0]
end

@testset "pruning a lone root removes it and compaction remaps from `start`" begin
    tree = LineageTree()
    r = P(tree, CE.NOPARENT, 1, 0.0, 1, 0, 0)
    CE._prune!(tree, r)
    @test CE._isdead(tree, r)
    @test isempty(CE._compact!(tree, Int32(1)) |> x -> filter(!=(0), x)) && length(tree) == 0

    tree = tree_fixture()
    CE._prune!(tree, Int32(3))               # root's right child
    remap = CE._compact!(tree, Int32(3))     # rows 1, 2 untouched
    @test tree.id == [1, 2, 4, 5]
    @test remap == [0, 3, 4]
    @test tree.left[2] == 3 && tree.right[2] == 4 && tree.right[1] == 0
    @test [l.id for l in Leaves(CellNode(tree, 1))] == [4, 5]
end

@testset "handles stay valid across compaction; removed cells throw" begin
    tree = tree_fixture()
    h5 = CellNode(tree, 5)
    h3 = CellNode(tree, 3)
    CE._prune!(tree, Int32(3))
    CE._compact!(tree, Int32(1))
    @test isalive(h5) && h5.id == 5 && h5.total == 9 && h5.clone == 3
    @test getfield(h5, :idx) == 5                     # stale hint, re-resolved by id
    @test h5.parent.id == 2
    @test !isalive(h3)
    @test_throws ArgumentError h3.clone
    @test !isalive(CellNode(tree, 1))                 # an internal node is not alive
end

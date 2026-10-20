# splitmix64 finaliser. Used instead of `Base.hash`, which is not guaranteed stable across
# Julia versions, so that a base seed reproduces the same draws everywhere. Changing it
# invalidates stored draws. Same function as `_splitmix64` in NonMarkovEvolution.jl.
function _splitmix64(x::UInt64)
    x += 0x9e3779b97f4a7c15
    x = (x ⊻ (x >> 30)) * 0xbf58476d1ce4e5b9
    x = (x ⊻ (x >> 27)) * 0x94d049bb133111eb
    return x ⊻ (x >> 31)
end

"""
    derive_seed(base::Integer, parts::Integer...) -> UInt64

Derive a reproducible `UInt64` seed from `base` and any number of integer `parts` by
chaining splitmix64, `h ← splitmix64(h ⊻ part)`, starting from `splitmix64(base)`. Stable
across Julia versions; distinct `parts` give independent-looking seeds.
"""
function derive_seed(base::Integer, parts::Integer...)
    h = _splitmix64(UInt64(base % UInt64))
    for p in parts
        h = _splitmix64(h ⊻ UInt64(p % UInt64))
    end
    return h
end

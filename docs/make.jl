using Documenter
using ClonalEvolutionTx

DocMeta.setdocmeta!(ClonalEvolutionTx, :DocTestSetup,
                    :(using ClonalEvolutionTx); recursive = true)

# Edit-on-GitHub links need a git commit and a remote; skip them outside CI so that a local
# build also works in a checkout without commits.
ci = get(ENV, "CI", "false") == "true"
local_only = ci ? (;) : (; remotes = nothing)

makedocs(;
    sitename = "ClonalEvolutionTx.jl",
    modules = [ClonalEvolutionTx],
    authors = "Alexander Stein",
    pages = [
        "Home" => "index.md",
        "Concepts" => "concepts.md",
        "The simulation block" => "blocks.md",
        "Mutations" => "mutations.md",
        "Treatment and schedules" => "treatment.md",
        "Measurements and output" => "measurements.md",
        "Tree statistics" => "statistics.md",
        "Sampling" => "sampling.md",
        "Performance" => "performance.md",
        "Limitations and open questions" => "limitations.md",
        "API reference" => "api.md",
    ],
    doctest = true,
    checkdocs = :exports,
    format = Documenter.HTML(
        prettyurls = ci,
        repolink = "https://github.com/alexsteininfo/ClonalEvolutionTx.jl",
        size_threshold_warn = 300 * 1024,    # the API page lists internals too
        size_threshold = 600 * 1024,
    ),
    local_only...,
)

deploydocs(
    repo = "github.com/alexsteininfo/ClonalEvolutionTx.jl.git",
    devbranch = "main",
)

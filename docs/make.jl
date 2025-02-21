using Documenter
using GraphPartitioning  # Make sure it uses the locally developed version

# println(nameof(GraphPartitioning))  # Should print GraphPartitioning
# println("Defined functions: ", names(GraphPartitioning, all=true))

makedocs(
    sitename = "GraphPartitioning.jl",
    modules  = [GraphPartitioning],
    format   = Documenter.HTML(
        prettyurls = false,  # Ensure local testing works
        collapselevel = 1    # Controls sidebar depth (1 = show subsections)
    ),
    repo     = "https://github.com/lechekhabm/GraphPartitioning.jl",
    pages    = [
        "Home" => "index.md",
        "Documentation"  => "documentation.md"
    ],
    warnonly = true,  # Prevent build failure due to missing docs
)


deploydocs(
    repo="github.com/lechekhabm/GraphPartitioning.jl.git",
    branch="gh-pages",
    devbranch="docs",
    push_preview=true
)
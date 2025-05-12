using Documenter
using GraphLab  # Make sure it uses the locally developed version

# println("Defined functions: ", names(GraphPartitioning, all=true))

makedocs(
    sitename = "GraphLab.jl",
    modules  = [GraphLab],
    format   = Documenter.HTML(
        repolink = "https://github.com/lechekhabm/GraphLab.jl",
        collapselevel = 1,
        prettyurls = true,
    ),
    repo     = "https://github.com/lechekhabm/GraphLab.jl",
    pages    = [
        "Home" => "index.md",
        "Usage Guide" => "usage.md",
        "API Reference" => "api.md",
        "Developers API Reference" => "dev_api.md",
    ],
    warnonly = true,  # Prevent build failure due to missing docs
)


deploydocs(;
    repo="github.com/lechekhabm/GraphLab.jl.git",
    branch="gh-pages",
    devbranch="main",
    versions = ["stable", "v#.#.#", "dev"],  # enable stable badge support
    forcepush=true,  # Ensure it force-pushes
    deploy_config=Documenter.GitHubActions()  # Adjust if your default branch is different
)
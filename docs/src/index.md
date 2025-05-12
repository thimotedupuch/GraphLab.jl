# GraphLab.jl Documentation

`GraphLab.jl` is a toolbox for graph partitioning algorithms in Julia. It provides a framework to compare, benchmark, and analyze various graph partitioning techniques. The framework supports both non-recursive and recursive methods with geometric information, like coordinate and inertial bisection, as well as methods without geometric information, such as spectral bisection. It makes it easy to evaluate and compare the outputs of different methods, helping users understand their trade-offs and applications.

## Getting Started

### Installation
To install the package from GitHub and add it to your working environment, follow these steps:

1. Add the `GraphLab.jl` to you project using the Julia command:
   ```julia
   using Pkg
   Pkg.add(url="https://github.com/lechekhabm/GraphLab.jl")
   ```
   **Note:** If you have an SSH key set up with GitHub, you can install the package using SSH:
   ```julia
   Pkg.add(url="https://<YOUR_GITHUB_TOKEN>@github.com/lechekhabm/GraphLab.jl.git")
   ```
   Replace `<YOUR_GITHUB_TOKEN>` with a valid token that has repository read access.

3. Basic Example:
   ```julia
   using GraphLab
   A, coords = GraphLab.build_adjacency("network")
   p = GraphLab.part_spectral(A)
   GraphLab.draw_graph(A, coords, p, file_name="test.png")
   ```
   
  For further details on the package and its functions, see the [Usage Guide](@ref).

### Prerequisites and Dependencies

The package will automatically install the following dependencies: **Arpack**, **CairoMakie**, **Colors**, **Graphs**, **LinearAlgebra**, **Metis**, **SparseArrays**, **Statistics**, **AMD**, **GraphsMatching**, **JuMP**, and **Cbc**.

For additional usability, you may need the following packages: **DelimitedFiles**, **MAT**, **Plots**, and **PrettyTables**.

Refer to the source code and `Project.toml` files within the `examples/` directory for further details about dependencies and configurations. Additional explanations are provided later in this document.

## Future Work
* Support for higher-dimensional graphs and hypergraphs
* Automated detection of external software
* Performance optimizations
* Support for additional partitioning methods
* Enhanced documentation and tutorials
* Interactive visualization

## Citation

Please cite [XXXX](https://epubs.siam.org/doi/10.1137/21M1392231) in your publications if it helps your research:
```
@article{doi:10.1137/21M1392231,
 author = {Name},
 title = {Tile},
 journal = X},
 volume = {44},
 number = {3},
 pages = {C210-C236},
 year = {2022},
 doi = {10.1137/21M1392231}
}
```
See [https://arxiv.org/pdf/XXXX.pdf](https://arxiv.org/pdf/2202.06555.pdf) for an archived version of the article. 

### Authors
* [NAME NAME](https://www.linkedin.com/in/mlechekhab/) (Institute, University)

* [Malik Lechekhab](https://www.linkedin.com/in/mlechekhab/) (Institute of Computing, Università della Svizzera italiana)

* [Aryan Eftekhari](https://scholar.google.com/citations?user=GiugKBsAAAAJ&hl=en) (Institute of Computing, Università della Svizzera italiana)

### Acknowledgments  
We would like to thank the authors of [SGtSNEpi.jl](https://github.com/fcdimitr/SGtSNEpi.jl) for their work, which served as a reference for the graph plotting function in this project.
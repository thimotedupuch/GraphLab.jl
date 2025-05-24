# GraphLab.jl Documentation

`GraphLab.jl` is a Julia package designed to make graph partitioning approachable, interactive, and educational.
Rather than focusing purely on production-ready pipelines, `GraphLab.jl` provides a hands-on framework for exploring, experimenting, and teaching graph partitioning techniques.

It offers an accessible collection of algorithms, including:
* Coordinate, inertial, spectral, and random sphere bisection.
* Space-filling curve partitioning.
* Recursive partitioning and nested dissection.

Beyond algorithms, `GraphLab.jl` includes tools for:
* Generating and manipulating adjacency matrices.
* Evaluating partition quality.
* Benchmarking algorithms.
* Visualizing partitioned graphs.

`GraphLab.jl` also integrates with external partitioners, enabling easy method comparisons in a unified Julia environment.

`GraphLab.jl` aims to support the Julia graph community's educational and exploratory needs, helping learners and researchers engage with graph theory in an interactive and practical way.

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
   A, coords = GraphLab.grid_graph(10, 50, π/3)
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

If you use **GraphLab.jl** in your research or teaching, please consider citing our upcoming paper:
```
@article{GraphLab2025,
author = {Malik Lechekhab and collaborators},
title = {GraphLab.jl: A Julia Framework for Exploring Graph Partitioning},
journal = {In preparation},
year = {2025},
note = {Preprint will be made available at arXiv},
}
```
We will update the citation once the paper is published.

## Authors
* [Malik Lechekhab](https://www.linkedin.com/in/mlechekhab/) (Institute of Computing, Università della Svizzera italiana)
* [Dimosthenis Pasadakis](https://search.usi.ch/en/people/bfe7763cea5221d043f905ad414e1a8d/pasadakis-dimosthenis) (Institute of Computing, Università della Svizzera italiana)
* [Aryan Eftekhari](https://scholar.google.com/citations?user=GiugKBsAAAAJ&hl=en) (Institute of Computing, Università della Svizzera italiana)
* [Roger Käppeli](https://math.ethz.ch/research/applied-mathematics-numerical-analysis-scientific-computing/roger-kaeppeli.html) (Department of Mathematics, ETH Zürich)
* [Olaf Schenk](https://search.usi.ch/en/people/9a52a2fdb8d3d26ec16fb1569b590909/schenk-olaf) (Institute of Computing, Università della Svizzera italiana)

## Acknowledgments  
- Jeroen Baert, Libmorton: C++ Morton Encoding/Decoding Library, 2018. [Link](https://github.com/Forceflow/libmorton).
- Aparna Sasidharan, Aparna, and Snir, Space-filling Curves for Partitioning Adaptively Refined Meshes.
- Aparna Sasidharan, John R. Gilbert, Shang-Hua Teng, Yingzhou Li, A General Space-filling Curve Algorithm for Partitioning 2D Meshes, In Proceedings of the 2015 IEEE 17th International Conference on High Performance Computing and Communications (HPCC), 2015. DOI: 10.1109/HPCC-CSS-ICESS.2015.192.
- Simon Byrne, GilbertCurves.jl, [Link](https://github.com/CliMA/GilbertCurves.jl).
- Jakub Červený, gilbert, [Link](https://github.com/jakubcerveny/gilbert).
- Yingzhou Li, meshpart Toolbox, [Link](https://github.com/YingzhouLi/meshpart).
We would also like to thank the authors of [SGtSNEpi.jl](https://github.com/fcdimitr/SGtSNEpi.jl) for their work, which served as a reference for the graph plotting function in this project.
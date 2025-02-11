# GraphPartitioning.jl

This is the code repository for ``GraphPartitioning.jl``, a toolbox for graph partitioning algorithms in Julia. It provides a framework to compare, benchmark, and analyze various graph partitioning techniques. The framework supports both non-recursive and recursive methods with geometric information, like coordinate and inertial bisection, as well as methods without geometric information, such as spectral bisection. It makes it easy to evaluate and compare the outputs of different methods, helping users understand their trade-offs and applications.

Feel free to contribute or report issues on the GitHub repository for **GraphPartitioning.jl**. Happy partitioning!

## Getting Started

### Installation
To install the package from GitHub and add it to your working environment, follow these steps:

1. Add the ``GraphPartitioning.jl`` to you project using the Julia command:
   ```julia
   using Pkg
   Pkg.add(url="https://github.com/lechekhabm/GraphPartitioning.jl")
   ```
   **Note:** If you have an SSH key set up with GitHub, you can install the package using SSH:
   ```julia
   Pkg.add(url="https://<YOUR_GITHUB_TOKEN>@github.com/lechekhabm/GraphPartitioning.jl.git")
   ```
   Replace `<YOUR_GITHUB_TOKEN>` with a valid token that has repository read access.

3. Basic Example:
   ```julia
   using GraphPartitioning
   A, coords = GraphPartitioning.build_adjacency("network")
   p = GraphPartitioning.part_spectral(A)
   GraphPartitioning.draw_graph(A, coords, ,p file_name="test.png")
   ```
   
  For further details on the package and its functions, see [Paper](www.google.com).

### Prerequisites and Dependencies

The package will automatically install the following dependencies: **Arpack**, **CairoMakie**, **Colors**, **Graphs**, **LinearAlgebra**, **Metis**, **SparseArrays**, and **Statistics**.

For additional usability, you may need the following packages: **DelimitedFiles**, **MAT**, **Plots**, and **PrettyTables**.

Refer to the source code and `Project.toml` files within the `examples/` directory for further details about dependencies and configurations. Additional explanations are provided later in this document.

## Running Examples

The package includes example scripts in the `GraphPartitioning/examples/` directory. These scripts demonstrate how to use the package to benchmark and compare different graph partitioning methods.

To run these examples you need to clone the repo using the command:
   ```bash
   git clone https://github.com/lechekhabm/GraphPartitioning.jl
   ```

### Example 1: Comparing Partitioning Methods

**Script:** `GraphPartitioning.jl/examples/ex1.jl`

This example compares different graph partitioning methods, including:

- **Coordinate Bisection**
- **Inertial Bisection**
- **Spectral Bisection**
- **Bisection using METIS**

The script evaluates these methods on a series of different graphs, providing insights into their performance and effectiveness.

<!-- 2x2 Grid for Airfoil1 Examples -->
<table style="width:100%; border-collapse: collapse;">
  <tr>
    <!-- Coordinate Bisection -->
    <td style="width:50%; text-align: center; padding: 0;">
      <img src="examples/ex1_airfoil1_coordinate.png" alt="Coordinate Bisection" style="width:100%; margin:0; padding:0;">
      <p style="margin: 0; padding: 0;"><em>Coordinate Bisection</em></p>
    </td>
    <!-- Inertial Bisection -->
    <td style="width:50%; text-align: center; padding: 0;">
      <img src="examples/ex1_airfoil1_inertial.png" alt="Inertial Bisection" style="width:100%; margin:0; padding:0;">
      <p style="margin: 0; padding: 0;"><em>Inertial Bisection</em></p>
    </td>
  </tr>
  <tr>
    <!-- Spectral Bisection -->
    <td style="width:50%; text-align: center; padding: 0;">
      <img src="examples/ex1_airfoil1_spectral.png" alt="Spectral Bisection" style="width:100%; margin:0; padding:0;">
      <p style="margin: 0; padding: 0;"><em>Spectral Bisection</em></p>
    </td>
    <!-- METIS Bisection -->
    <td style="width:50%; text-align: center; padding: 0;">
      <img src="examples/ex1_airfoil1_metis.png" alt="METIS Bisection" style="width:100%; margin:0; padding:0;">
      <p style="margin: 0; padding: 0;"><em>METIS</em></p>
    </td>
  </tr>
</table>

#### How to Run Example 1

1. Navigate to the `examples` directory:

   ```bash
   cd GraphPartitioning.jl/examples
   ```

2. Run the example script in the Julia or directly from the terminal:

   ```bash
   julia ex1.jl
   ```

### Example 2: Recursive Bisection

**Script:** `GraphPartitioning.jl/examples/ex2.jl`

This example demonstrates recursive bisection using different methods, including:

- **Coordinate Bisection**
- **Inertial Bisection**
- **Spectral Bisection**
- **Recursive Bisection using METIS K-way**

<table style="width:100%; border-collapse: collapse;">
  <tr>
    <!-- Recursive Coordinate Bisection -->
    <td style="width:50%; text-align: center; padding: 10px;">
      <img src="examples/ex2_Swiss_graph_coordinate.png" alt="Recursive Coordinate Bisection" style="width:100%;">
      <p style="font-size: 12px;"><em>Recursive Coordinate Bisection</em></p>
    </td>
    <!-- Recursive Inertial Bisection -->
    <td style="width:50%; text-align: center; padding: 10px;">
      <img src="examples/ex2_Swiss_graph_inertial.png" alt="Recursive Inertial Bisection" style="width:100%;">
      <p style="font-size: 12px;"><em>Recursive Inertial Bisection</em></p>
    </td>
  </tr>
  <tr>
    <!-- Recursive Spectral Bisection -->
    <td style="width:50%; text-align: center; padding: 10px;">
      <img src="examples/ex2_Swiss_graph_spectral.png" alt="Recursive Spectral Bisection" style="width:100%;">
      <p style="font-size: 12px;"><em>Recursive Spectral Bisection</em></p>
    </td>
    <!-- Recursive METIS Bisection -->
    <td style="width:50%; text-align: center; padding: 10px;">
      <img src="examples/ex2_Swiss_graph_metis_rec.png" alt="Recursive METIS Bisection" style="width:100%;">
      <p style="font-size: 12px;"><em>Recursive METIS Bisection (K-way is also available)</em></p>
    </td>
  </tr>
</table>

#### How to Run Example 2

1. Navigate to the `examples` directory:

   ```bash
   cd GraphPartitioning.jl/examples
   ```

2. Run the example script in the Julia or directly from the terminal:

   ```bash
   julia ex2.jl
   ```

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
See [https://arxiv.org/pdf/2202.06555.pdf](https://arxiv.org/pdf/2202.06555.pdf) for an archived version of the article. 

### Authors
* [Name Name](https://google.com) (Department, University of X)

* [Aryan Eftekhari](https://scholar.google.com/citations?user=GiugKBsAAAAJ&hl=en) (Institute of Computing, Università della Svizzera italiana)

### Acknowledgments  
We would like to thank the authors of [SGtSNEpi.jl](https://github.com/fcdimitr/SGtSNEpi.jl) for their work, which served as a reference for the graph plotting function in this project.

# Usage Guide

This guide shows how to use `GraphLab.jl` for various tasks.

## Running Examples

The package includes example scripts in the `GraphLab/examples/` directory. These scripts demonstrate how to use the package to benchmark and compare different graph partitioning methods.

To run these examples you need to clone the repo using the command:
```shell
git clone https://github.com/lechekhabm/GraphLab.jl
```

### Quick Start Example
```julia
using GraphLab

A, coords = GraphLab.build_adjacency("network")
p = GraphLab.part_spectral(A)
GraphLab.draw_graph(A, coords, p, file_name="test.png")
```

A simple example demonstrating the usage of `GraphLab.jl` is available in this **Google Colab notebook**: [Open in Colab](https://colab.research.google.com/drive/1cxkk3HV11BRfoOp27GIRVuL9mdzzZIjP?usp=sharing).

### Example 1: Comparing Partitioning Methods

**Script:** `GraphLab.jl/examples/ex1.jl`

This example compares different graph partitioning methods, including:

- **Coordinate Bisection**
- **Inertial Bisection**
- **Spectral Bisection**
- **Bisection using METIS**

The script evaluates these methods on a series of different graphs, providing insights into their performance and effectiveness.
```@raw html
<!-- 2x2 Grid for Airfoil1 Examples -->
<table style="width:100%; border-collapse: collapse;">
  <tr>
    <!-- Coordinate Bisection -->
    <td style="width:50%; text-align: center; padding: 0;">
      <img src="https://github.com/lechekhabm/GraphLab.jl/blob/main/examples/ex1_airfoil1_coordinate.png?raw=true" alt="Coordinate Bisection" style="width:100%; margin:0; padding:0;">
      <p style="margin: 0; padding: 0;"><em>Coordinate Bisection</em></p>
    </td>
    <!-- Inertial Bisection -->
    <td style="width:50%; text-align: center; padding: 0;">
      <img src="https://github.com/lechekhabm/GraphLab.jl/blob/main/examples/ex1_airfoil1_inertial.png?raw=true" alt="Inertial Bisection" style="width:100%; margin:0; padding:0;">
      <p style="margin: 0; padding: 0;"><em>Inertial Bisection</em></p>
    </td>
  </tr>
  <tr>
    <!-- Spectral Bisection -->
    <td style="width:50%; text-align: center; padding: 0;">
      <img src="https://github.com/lechekhabm/GraphLab.jl/blob/main/examples/ex1_airfoil1_spectral.png?raw=true" alt="Spectral Bisection" style="width:100%; margin:0; padding:0;">
      <p style="margin: 0; padding: 0;"><em>Spectral Bisection</em></p>
    </td>
    <!-- METIS Bisection -->
    <td style="width:50%; text-align: center; padding: 0;">
      <img src="https://github.com/lechekhabm/GraphLab.jl/blob/main/examples/ex1_airfoil1_metis.png?raw=true" alt="METIS Bisection" style="width:100%; margin:0; padding:0;">
      <p style="margin: 0; padding: 0;"><em>METIS</em></p>
    </td>
  </tr>
</table>
```

#### How to Run Example 1

1. Navigate to the `examples` directory:

   ```bash
   cd GraphLab.jl/examples
   ```

2. Run the example script in the Julia or directly from the terminal:

   ```bash
   julia ex1.jl
   ```

### Example 2: Recursive Bisection

**Script:** `GraphLab.jl/examples/ex2.jl`

This example demonstrates recursive bisection using different methods, including:

- **Coordinate Bisection**
- **Inertial Bisection**
- **Spectral Bisection**
- **Recursive Bisection using METIS K-way**

```@raw html
<table style="width:100%; border-collapse: collapse;">
  <tr>
    <!-- Recursive Coordinate Bisection -->
    <td style="width:50%; text-align: center; padding: 10px;">
      <img src="https://github.com/lechekhabm/GraphLab.jl/blob/main/examples/ex2_Swiss_graph_coordinate.png?raw=true" alt="Recursive Coordinate Bisection" style="width:100%;">
      <p style="font-size: 12px;"><em>Recursive Coordinate Bisection</em></p>
    </td>
    <!-- Recursive Inertial Bisection -->
    <td style="width:50%; text-align: center; padding: 10px;">
      <img src="https://github.com/lechekhabm/GraphLab.jl/blob/main/examples/ex2_Swiss_graph_inertial.png?raw=true" alt="Recursive Inertial Bisection" style="width:100%;">
      <p style="font-size: 12px;"><em>Recursive Inertial Bisection</em></p>
    </td>
  </tr>
  <tr>
    <!-- Recursive Spectral Bisection -->
    <td style="width:50%; text-align: center; padding: 10px;">
      <img src="https://github.com/lechekhabm/GraphLab.jl/blob/main/examples/ex2_Swiss_graph_spectral.png?raw=true" alt="Recursive Spectral Bisection" style="width:100%;">
      <p style="font-size: 12px;"><em>Recursive Spectral Bisection</em></p>
    </td>
    <!-- Recursive METIS Bisection -->
    <td style="width:50%; text-align: center; padding: 10px;">
      <img src="https://github.com/lechekhabm/GraphLab.jl/blob/main/examples/ex2_Swiss_graph_metis_rec.png?raw=true" alt="Recursive METIS Bisection" style="width:100%;">
      <p style="font-size: 12px;"><em>Recursive METIS Bisection (K-way is also available)</em></p>
    </td>
  </tr>
</table>
```

#### How to Run Example 2

1. Navigate to the `examples` directory:

   ```bash
   cd GraphLab.jl/examples
   ```

2. Run the example script in the Julia or directly from the terminal:

   ```bash
   julia ex2.jl
   ```
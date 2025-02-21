# Documentation
## Graph Construction
```@docs
GraphPartitioning.build_adjacency
```
## Native Partitioning Methods
```@docs
GraphPartitioning.part_coordinate
GraphPartitioning.part_inertial
GraphPartitioning.part_spectral
GraphPartitioning.part_metis
```
## External Partitioning Methods

The following functions interface with external graph partitioning software and require that the corresponding tools be installed on your system. Additionally, you must ensure that the environment variables for these executables are correctly set before running the functions.

---
### Graclus (`run_graclus`) 
Graclus is a spectral-based graph clustering and partitioning software.  

- **Installation**: You need to manually **download and compile Graclus** from its [official source](https://www.cs.utexas.edu/~dml/Software/graclus.html).  
- **Environment Variable**: Set `GRACLUS_PATH` to point to the compiled Graclus executable.  

```bash
export GRACLUS_PATH="/path/to/graclus"
```
```@docs
GraphPartitioning.run_graclus
```
#### Citation
* Weighted Graph Cuts without Eigenvectors: A Multilevel Approach, I. Dhillon, Y. Guan, and B, Kulis, IEEE Transactions on Pattern Analysis and Machine Intelligence (PAMI), vol. 29:11, pages 1944-1957, November 2007.
* A Fast Kernel-based Multilevel Algorithm for Graph Clustering, I. Dhillon, Y. Guan, and B, Kulis, Proceedings of The 11th ACM SIGKDD, Chicago, IL, August 21-24, 2005.
* Kernel k-means, Spectral Clustering and Normalized Cuts, I. Dhillon, Y. Guan, and B. Kulis, Proceedings of The 10th ACM SIGKDD, Seattle, WA, August 22-25, 2004.
---

### KaHIP `(run_kahip)`
KaHIP (Karlsruhe High-Quality Partitioning) is an efficient graph partitioning tool that supports various heuristics and optimizations.

- **Installation**: You need to **download and compile KaHIP** from its [official source](https://kahip.github.io/).
- **Environment Variable**: Set `KAHIP_PATH` to the KaHIP binary.

```bash
export KAHIP_PATH="/path/to/kahip"
```
```@docs
GraphPartitioning.run_kahip
```
#### Citation
* Peter Sanders and Christian Schulz. Engineering Multilevel Graph Partitioning Algorithms. In Proceedings of the 19th European Symposium on Algorithms (ESA'11), volume 6942 of LNCS, pages 469--480. Springer, 2011.
---
### Multi-Level Diffusion Clustering (MDC) `(run_mdc)`
MDC is a novel diffusion-based graph partitioning method.

- **Installation**: You must **download, compile, and install MDC** following the instructions from [official source](https://kahip.github.io/).
- **Environment Variable**: Set `MDC_PATH` to point to the MDC executable. It also requires setting `LD_LIBRARY_PATH` for proper execution.
```bash
export MDC_PATH="/path/to/mdc"
export LD_LIBRARY_PATH="/path/to/stag/stag-1.3.0/build_dir/stag_lib/"
```
```@docs
GraphPartitioning.run_mdc
```
#### Citation
* Lechekhab M., Pasadakis D., Schenk O. (forthcoming) Multilevel Diffusion Based Spectral Graph Clustering. 28th Annual IEEE High Performance Extreme Computing Virtual Conference. IEEE. 28th Annual IEEE High Performance Extreme Computing Virtual Conference. Virtual Conference. September 23-27, 2024
---
## Recursive Bisection
```@docs
GraphPartitioning.recursive_bisection
```
## Visualization
```@docs
GraphPartitioning.draw_graph
```

## Utility Function
```@docs
GraphPartitioning.count_edge_cut
```

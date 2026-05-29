# persistent-sheaf

Persistent sheaf cohomology, cellular sheaf Laplacians, and multi-modal data fusion via topological data analysis.

## Usage

```rust
use persistent_sheaf::{SimplicialComplex, PersistenceDiagram, CellularSheaf, Filtration};

// Build a Vietoris-Rips complex from distances
let distances = vec![vec![0.0, 1.0, 2.0], vec![1.0, 0.0, 1.0], vec![2.0, 1.0, 0.0]];
let complex = SimplicialComplex::vietoris_rips(&distances, 1.5);

// Compute persistence
let filtration = Filtration::from_distance_matrix(&distances, 10);
let diagram = filtration.compute_persistence();

// Cellular sheaf with Laplacian
let sheaf = CellularSheaf::constant(complex, 2);
let laplacian = SheafLaplacian::from_sheaf(&sheaf);
```

## Features

- **Simplicial complexes**: Vertices, edges, triangles with Euler characteristic and Betti numbers
- **Vietoris-Rips construction** from distance matrices
- **Persistence diagrams**: Birth-death pairs, bottleneck distance, Betti curves
- **Cellular sheaves**: Constant and weighted with cohomology computation
- **Sheaf Laplacian**: Generalizes graph Laplacian with sheaf-theoretic information
- **Filtration builder**: From distance matrices with incremental persistence

## Tests

28 tests, all passing. `cargo test` to run.

## License

MIT

Part of the [SuperInstance OpenConstruct](https://github.com/SuperInstance/OpenConstruct) ecosystem.

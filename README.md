# persistent-sheaf

Topological data analysis via persistent sheaf cohomology. Point this at any dataset and it tells you the shape of your data.

## The Shape of Data

You have points. Clusters, holes, tunnels, voids — the *topology* of your data is real information. Persistent homology finds that shape. Sheaf theory tells you how different measurements of the same system relate. This crate does both.

```rust
use persistent_sheaf::*;
```

---

## 1. Build a Complex from Real Data

Five temperature sensors spread across a building. Two clusters: the server room runs hot, the offices run cool. One sensor in the hallway between them.

```rust
use persistent_sheaf::SimplicialComplex;

// Sensor positions: [x, y] — 2 in the server room, 2 in the office wing, 1 in between
let distances = vec![
    // S0   S1   S2   S3   S4
    vec![0.0, 0.8, 4.0, 4.2, 2.0], // S0: server room
    vec![0.8, 0.0, 3.5, 3.8, 1.8], // S1: server room
    vec![4.0, 3.5, 0.0, 0.9, 1.9], // S2: office wing
    vec![4.2, 3.8, 0.9, 0.0, 2.1], // S3: office wing
    vec![2.0, 1.8, 1.9, 2.1, 0.0], // S4: hallway
];

let complex = SimplicialComplex::vietoris_rips(&distances, 2.0);

println!("Vertices: {}", complex.num_simplices(0)); // 5
println!("Edges:    {}", complex.num_simplices(1)); // edges where dist ≤ 2.0
println!("Triangles: {}", complex.num_simplices(2));
println!("Euler characteristic: {}", complex.euler_characteristic());

let betti = complex.betti_numbers();
println!("β₀ = {} (connected components — your clusters)", betti[0]);
println!("β₁ = {} (1-dimensional holes/loops)", betti[1]);
// With epsilon=2.0: S0-S1, S4 connected to both groups
// β₀ = 1 (one component), β₁ depends on whether triangles fill the loop
```

At `epsilon = 1.0`, you get two clusters (β₀=2). At `epsilon = 2.0`, the hallway connects them (β₀=1). That transition IS the information.

---

## 2. Watch Topology Appear and Disappear

A filtration builds the complex at every scale. Persistence diagrams record which features survive.

```rust
use persistent_sheaf::Filtration;
use persistent_sheaf::SimplicialComplex;

// The same sensor network
let distances = vec![
    vec![0.0, 0.8, 4.0, 4.2, 2.0],
    vec![0.8, 0.0, 3.5, 3.8, 1.8],
    vec![4.0, 3.5, 0.0, 0.9, 2.1],
    vec![4.2, 3.8, 0.9, 0.0, 2.3],
    vec![2.0, 1.8, 2.1, 2.3, 0.0],
];

// Sweep epsilon from 0 to max distance in 20 steps
let filtration = Filtration::from_distance_matrix(&distances, 20);
let diagram = filtration.compute_persistence();

println!("Found {} persistence pairs", diagram.len());

// Read the diagram
for pair in &diagram.pairs {
    let label = match pair.dimension {
        0 => "CLUSTER",
        1 => "HOLE",
        _ => "FEATURE",
    };
    let persist = pair.persistence();
    if pair.is_essential() {
        println!("[{}] dim={} birth={:.2} death=∞  — essential (never dies)", label, pair.dimension, pair.birth);
    } else {
        println!("[{}] dim={} birth={:.2} death={:.2}  persistence={:.2}",
            label, pair.dimension, pair.birth, pair.death, persist);
    }
}

// Long bars = real structure. Short bars = noise.
// The cluster merging at epsilon≈2.0? That's a real bar.
// A tiny 1-dimensional hole at epsilon=1.5 that dies at epsilon=1.7? Noise.
```

**How to read the diagram:**
- Long bar at dimension 0 → a real cluster in your data
- Long bar at dimension 1 → a real loop/tunnel
- Short bars → noise, not structure
- Essential features (death=∞) → survive at all scales

---

## 3. Find the Most Important Feature

```rust
use persistent_sheaf::PersistenceDiagram;

let mut diagram = PersistenceDiagram::new();
// From our sensor network analysis
diagram.add(0.0, 0.8, 0);   // cluster 1: appears immediately, merges at 0.8
diagram.add(0.0, 1.8, 0);   // cluster 2: appears, merges when hallway connects
diagram.add(0.0, f64::INFINITY, 0); // the final unified cluster (essential)
diagram.add(1.5, 2.0, 1);   // a brief 1D hole (noise)

// What's the most persistent feature?
if let Some(mp) = diagram.most_persistent() {
    println!("Most persistent: dim={}, persists for {:.2}",
        mp.dimension, mp.persistence());
}

// How much total structure? (higher power = emphasizes long bars more)
let total_p1 = diagram.total_persistence(1.0);
let total_p2 = diagram.total_persistence(2.0);
println!("Total persistence (p=1): {:.2}", total_p1);
println!("Total persistence (p=2): {:.2} — squared penalizes short bars", total_p2);

// Filter to just the clusters
let dim0 = diagram.filter_dimension(0);
println!("{} cluster-level features", dim0.len());

// Compare two datasets: are they topologically similar?
let mut diagram2 = PersistenceDiagram::new();
diagram2.add(0.0, f64::INFINITY, 0);
let bottleneck = diagram.bottleneck_distance(&diagram2);
println!("Bottleneck distance to a single cluster: {:.2}", bottleneck);
// Small distance = similar topology. Large = different shape.
```

---

## 4. Betti Curves: How Topology Changes With Scale

```rust
use persistent_sheaf::PersistenceDiagram;

let mut diagram = PersistenceDiagram::new();
// Three clusters merge at different scales
diagram.add(0.0, 1.0, 0);   // first merge
diagram.add(0.0, 2.0, 0);   // second merge
diagram.add(0.0, f64::INFINITY, 0); // survives forever
diagram.add(1.5, 2.5, 1);   // a 1D hole that appears then fills

let thresholds: Vec<f64> = (0..50).map(|i| i as f64 * 0.1).collect();
let betti_curve = diagram.betti_curve(&thresholds);

println!("Betti curve (number of alive features vs. scale):");
for (i, &t) in thresholds.iter().enumerate() {
    if betti_curve[i] > 0 {
        println!("  ε={:.1}: {} features alive", t, betti_curve[i]);
    }
}
// At ε=0.0: 3 features (all three clusters)
// At ε=0.5: 3 features (nothing merged yet)
// At ε=1.0: 2 features (first cluster merged)
// At ε=1.5: 2 + 1 hole = 3 features (hole appeared!)
// At ε=2.0: 1 + 1 = 2 features (second merge, hole still alive)
// At ε=2.5: 1 feature (only the essential one survives)
```

The Betti curve is a fingerprint. Same-shaped data produces the same curve.

---

## 5. Cellular Sheaves: When Different Sensors See Different Things

A sheaf assigns data to each cell and tracks how they must agree. The sheaf Laplacian tells you where they don't.

```rust
use persistent_sheaf::{SimplicialComplex, CellularSheaf, SheafLaplacian};

// Build the complex: 3 agents monitoring a system
let mut complex = SimplicialComplex::new();
complex.add_edge(0, 1); // agent 0 ↔ agent 1
complex.add_edge(1, 2); // agent 1 ↔ agent 2

// Constant sheaf: every agent has the same 2D state space,
// restriction maps are identity (perfect agreement required)
let sheaf = CellularSheaf::constant(complex.clone(), 2);

println!("Stalk dimension: {} (each agent has R² state)", sheaf.stalk_dimension);
println!("Global section dimension: {}", sheaf.global_section_dimension());
// Global sections = assignments where all agents agree = 2 (the full R²)
println!("H⁰ dimension: {}", sheaf.cohomology_dimension(0));
println!("H¹ dimension: {}", sheaf.cohomology_dimension(1));
// H⁰ > 0 means global agreement is possible
// H¹ > 0 means there's obstruction to extending local info globally

// Build the sheaf Laplacian
let laplacian = SheafLaplacian::from_sheaf(&sheaf);
println!("Laplacian dimension: {}×{}", laplacian.dimension, laplacian.dimension);
// 3 vertices × 2 stalk dim = 6×6 matrix

// The sheaf Laplacian generalizes the graph Laplacian.
// Its kernel = global sections (where all agents agree).
// Its eigenvalues = how "rigid" the agreement constraints are.
let largest_eig = laplacian.largest_eigenvalue(100);
println!("Largest eigenvalue: {:.4}", largest_eig);
// Large eigenvalue = strong disagreement penalty = rigid system
```

Now with *weighted* sheaf — different agents have different trust levels:

```rust
use persistent_sheaf::{SimplicialComplex, CellularSheaf, SheafLaplacian};

let mut complex = SimplicialComplex::new();
complex.add_edge(0, 1);
complex.add_edge(1, 2);
complex.add_edge(0, 2); // triangle — fully connected

// Weighted sheaf: edge 0→1 is a strong link, 1→2 is weak, 0→2 is medium
let sheaf = CellularSheaf::from_weights(complex, &[1.0, 0.3, 0.7]);

let laplacian = SheafLaplacian::from_sheaf(&sheaf);
println!("Weighted Laplacian: {}×{}", laplacian.dimension, laplacian.dimension);

// Compare with standard graph Laplacian
let graph_lap = SheafLaplacian::graph_laplacian(3, &[(0, 1), (1, 2), (0, 2)]);
let fiedler = graph_lap.fiedler_value();
println!("Fiedler value (algebraic connectivity): {:.4}", fiedler);
// Higher = more connected graph = faster information diffusion
```

---

## 6. Graph Laplacian: The Special Case

When all sheaf stalks are ℝ and all restriction maps are 1, the sheaf Laplacian IS the graph Laplacian.

```rust
use persistent_sheaf::SheafLaplacian;

// A path graph: 0 — 1 — 2 — 3 — 4
let edges = vec![(0, 1), (1, 2), (2, 3), (3, 4)];
let lap = SheafLaplacian::graph_laplacian(5, &edges);

// The Laplacian matrix
println!("Graph Laplacian (path on 5 vertices):");
for i in 0..5 {
    print!("  [");
    for j in 0..5 {
        print!("{:5.1}", lap.matrix[i][j]);
    }
    println!("]");
}
// [ 1.0 -1.0  0.0  0.0  0.0]
// [-1.0  2.0 -1.0  0.0  0.0]
// [ 0.0 -1.0  2.0 -1.0  0.0]
// [ 0.0  0.0 -1.0  2.0 -1.0]
// [ 0.0  0.0  0.0 -1.0  1.0]

// Multiply by a signal
let signal = vec![1.0, 0.0, 0.0, 0.0, 0.0]; // impulse at vertex 0
let smooth = lap.mul_vec(&signal);
println!("L * [1,0,0,0,0] = {:?}", smooth);
// [1, -1, 0, 0, 0] — the "roughness" of the signal

let eig = lap.largest_eigenvalue(100);
println!("Largest eigenvalue: {:.4}", eig);
// For a path graph this is ~3.62
```

---

## 7. Full Pipeline: From Raw Data to Topology Report

```rust
use persistent_sheaf::*;

// Agent performance metrics over 6 runs
// Each row = [latency_ms, throughput, error_rate]
let metrics: Vec<Vec<f64>> = vec![
    vec![12.0, 950.0, 0.01], // run 0: good
    vec![13.0, 940.0, 0.02], // run 1: good
    vec![11.0, 960.0, 0.01], // run 2: good
    vec![45.0, 500.0, 0.15], // run 3: degraded
    vec![50.0, 480.0, 0.18], // run 4: degraded
    vec![12.0, 945.0, 0.01], // run 5: good
];

// Compute distance matrix (Euclidean)
let n = metrics.len();
let distances: Vec<Vec<f64>> = (0..n)
    .map(|i| (0..n)
        .map(|j| metrics[i].iter().zip(&metrics[j])
            .map(|(a, b)| (a - b).powi(2)).sum::<f64>().sqrt())
        .collect())
    .collect();

// Step 1: Filtration — sweep across all scales
let filtration = Filtration::from_distance_matrix(&distances, 30);
println!("Filtration has {} steps", filtration.len());

// Step 2: Persistence — which features are real?
let diagram = filtration.compute_persistence();
println!("\n=== Topology Report ===");
for pair in &diagram.pairs {
    if pair.persistence() > 5.0 {
        println!("  ROBUST: dim={} persists {:.1} units (birth={:.1}, death={:.1})",
            pair.dimension, pair.persistence(), pair.birth,
            if pair.is_essential() { f64::INFINITY } else { pair.death });
    }
}

// Step 3: Betti curve — how many features at each scale?
let thresholds: Vec<f64> = (0..20).map(|i| i as f64 * 5.0).collect();
let curve = diagram.betti_curve(&thresholds);
let max_features = curve.iter().max().copied().unwrap_or(0);
println!("\nPeak feature count: {}", max_features);

// Step 4: Sheaf analysis — are the runs locally consistent?
let mut complex = SimplicialComplex::vietoris_rips(&distances, 20.0);
let sheaf = CellularSheaf::constant(complex, 1); // 1D stalk = scalar metric
let laplacian = SheafLaplacian::from_sheaf(&sheaf);
println!("\nSheaf Laplacian: {}×{}", laplacian.dimension, laplacian.dimension);
println!("Largest eigenvalue: {:.2}", laplacian.largest_eigenvalue(100));
// The eigenvalue spectrum tells you how "rigid" the cluster assignments are
```

Point this at any distance matrix. It finds the clusters, the holes, the noise. The Betti curve is your topology fingerprint. The sheaf Laplacian tells you how rigidly your data points are constrained to agree.

---

## API Reference

| Type | What it does |
|------|-------------|
| `SimplicialComplex` | Vertices, edges, triangles, tetrahedra. Betti numbers. Vietoris-Rips construction. |
| `Filtration` | Nested sequence of complexes built by sweeping epsilon. Computes persistence. |
| `PersistenceDiagram` | Birth-death pairs. Most persistent feature. Betti curves. Bottleneck distance. |
| `CellularSheaf` | Assigns vector spaces (stalks) to cells with restriction maps. Cohomology. |
| `SheafLaplacian` | The operator that generalizes graph Laplacian with sheaf data. Eigenvalues = rigidity. |

## The Pipeline

```
raw data → distance matrix → Vietoris-Rips complex → filtration → persistence diagram
                                                                    ↕
                                              Betti curves ← ← ← ← ←
                                                                    
complex + restriction maps → cellular sheaf → sheaf Laplacian → eigenvalue spectrum
```

Every step is a function. Every function has a clear input and output. Compose them.

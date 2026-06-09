//! # Persistent Sheaf Tutorial
//!
//! A progressive introduction to topological data analysis with persistent sheaf cohomology.
//!
//! Run with: `cargo run --example tutorial`

use persistent_sheaf::{
    CellularSheaf, Filtration, PersistenceDiagram, SheafLaplacian, SimplicialComplex,
};

fn main() {
    lesson_1_building_a_complex();
    lesson_2_euler_characteristic();
    lesson_3_vietoris_rips();
    lesson_4_persistence_diagrams();
    lesson_5_filtration_and_homology();
    lesson_6_betti_curve_and_bottleneck();
    lesson_7_cellular_sheaves();
    lesson_8_sheaf_laplacian();
}

// ── Lesson 1: Building a Simplicial Complex ─────────────────────────────
fn lesson_1_building_a_complex() {
    println!("═══════════════════════════════════════════");
    println!("Lesson 1: Building a Simplicial Complex");
    println!("═══════════════════════════════════════════");
    println!();

    // Start with an empty complex and add simplices one at a time.
    let mut cx = SimplicialComplex::new();
    cx.add_vertex(0);
    cx.add_vertex(1);
    cx.add_vertex(2);
    println!("After adding 3 vertices:");
    println!("  0-simplices (vertices): {}", cx.num_simplices(0));
    println!("  1-simplices (edges):    {}", cx.num_simplices(1));
    println!();

    // Adding an edge automatically ensures its vertices exist.
    cx.add_edge(0, 1);
    cx.add_edge(1, 2);
    cx.add_edge(0, 2);
    println!("After adding 3 edges:");
    println!("  0-simplices: {}", cx.num_simplices(0));
    println!("  1-simplices: {}", cx.num_simplices(1));
    println!();

    // A triangle creates its boundary edges and vertices.
    let mut tri = SimplicialComplex::new();
    tri.add_triangle(0, 1, 2);
    println!("After add_triangle(0,1,2):");
    println!("  vertices: {}, edges: {}, triangles: {}",
        tri.num_simplices(0), tri.num_simplices(1), tri.num_simplices(2));
    println!();
}

// ── Lesson 2: Euler Characteristic ──────────────────────────────────────
fn lesson_2_euler_characteristic() {
    println!("═══════════════════════════════════════════");
    println!("Lesson 2: Euler Characteristic and Betti Numbers");
    println!("═══════════════════════════════════════════");
    println!();

    // χ = V − E + F − T (Euler's formula).
    let mut triangle = SimplicialComplex::new();
    triangle.add_triangle(0, 1, 2);
    println!("Filled triangle (V=3, E=3, F=1):");
    println!("  χ = {} (should be 1)", triangle.euler_characteristic());
    println!();

    // Hollow triangle (just edges, no face).
    let mut hollow = SimplicialComplex::new();
    hollow.add_edge(0, 1);
    hollow.add_edge(1, 2);
    hollow.add_edge(0, 2);
    println!("Hollow triangle (V=3, E=3, F=0):");
    println!("  χ = {} (should be 0)", hollow.euler_characteristic());
    println!();

    // Betti numbers: β₀ = connected components, β₁ = independent cycles.
    let betti_hollow = hollow.betti_numbers();
    println!("Hollow triangle Betti numbers: β₀={}, β₁={}", betti_hollow[0], betti_hollow[1]);
    let betti_filled = triangle.betti_numbers();
    println!("Filled triangle Betti numbers: β₀={}, β₁={}", betti_filled[0], betti_filled[1]);
    println!();
}

// ── Lesson 3: Vietoris-Rips Complex ─────────────────────────────────────
fn lesson_3_vietoris_rips() {
    println!("═══════════════════════════════════════════");
    println!("Lesson 3: Vietoris-Rips Complex");
    println!("═══════════════════════════════════════════");
    println!();

    // Build a VR complex from a distance matrix at varying thresholds.
    let distances = vec![
        vec![0.0, 1.0, 2.0, 3.0],
        vec![1.0, 0.0, 1.0, 2.0],
        vec![2.0, 1.0, 0.0, 1.0],
        vec![3.0, 2.0, 1.0, 0.0],
    ];

    // ε = 0.5: all points isolated.
    let small = SimplicialComplex::vietoris_rips(&distances, 0.5);
    println!("ε = 0.5: vertices={}, edges={}", small.num_simplices(0), small.num_simplices(1));

    // ε = 1.0: adjacent points connected.
    let medium = SimplicialComplex::vietoris_rips(&distances, 1.0);
    println!("ε = 1.0: vertices={}, edges={}, triangles={}",
        medium.num_simplices(0), medium.num_simplices(1), medium.num_simplices(2));

    // ε = 3.0: everything connected.
    let large = SimplicialComplex::vietoris_rips(&distances, 3.0);
    println!("ε = 3.0: vertices={}, edges={}, triangles={}",
        large.num_simplices(0), large.num_simplices(1), large.num_simplices(2));
    println!();
}

// ── Lesson 4: Persistence Diagrams ──────────────────────────────────────
fn lesson_4_persistence_diagrams() {
    println!("═══════════════════════════════════════════");
    println!("Lesson 4: Persistence Diagrams");
    println!("═══════════════════════════════════════════");
    println!();

    // A persistence diagram records birth–death of topological features.
    let mut dg = PersistenceDiagram::new();
    dg.add(0.0, 1.0, 0);   // H₀ feature born at 0, dies at 1
    dg.add(0.0, 2.0, 0);   // longer-lived H₀
    dg.add(0.5, 1.5, 1);   // H₁ loop born at 0.5, dies at 1.5
    println!("Diagram has {} pairs", dg.len());
    println!();

    // Inspect individual pairs.
    println!("H₀ features:");
    for p in dg.filter_dimension(0) {
        println!("  birth={:.1}, death={:.1}, persistence={:.1}, essential={}",
            p.birth, p.death, p.persistence(), p.is_essential());
    }
    println!("H₁ features:");
    for p in dg.filter_dimension(1) {
        println!("  birth={:.1}, death={:.1}, midpoint={:.2}",
            p.birth, p.death, p.midpoint());
    }
    println!();

    // Most persistent feature.
    if let Some(mp) = dg.most_persistent() {
        println!("Most persistent: birth={:.1}, death={:.1}, dim={}", mp.birth, mp.death, mp.dimension);
    }
    println!();
}

// ── Lesson 5: Filtrations and Persistent Homology ───────────────────────
fn lesson_5_filtration_and_homology() {
    println!("═══════════════════════════════════════════");
    println!("Lesson 5: Filtrations and Persistent Homology");
    println!("═══════════════════════════════════════════");
    println!();

    // A filtration builds VR complexes at increasing scales.
    let distances = vec![
        vec![0.0, 1.0, 2.0],
        vec![1.0, 0.0, 1.0],
        vec![2.0, 1.0, 0.0],
    ];
    let filt = Filtration::from_distance_matrix(&distances, 5);
    println!("Filtration has {} steps", filt.len());
    for (threshold, cx) in &filt.complexes {
        println!("  ε={:.2}: vertices={}, edges={}, triangles={}",
            threshold, cx.num_simplices(0), cx.num_simplices(1), cx.num_simplices(2));
    }
    println!();

    // Compute persistent homology from the filtration.
    let dg = filt.compute_persistence();
    println!("Persistence diagram: {} pairs", dg.len());
    for p in &dg.pairs {
        let death = if p.death == f64::INFINITY { "∞".to_string() } else { format!("{:.2}", p.death) };
        println!("  H{}: birth={:.2} → death={}", p.dimension, p.birth, death);
    }
    println!();
}

// ── Lesson 6: Betti Curves and Bottleneck Distance ──────────────────────
fn lesson_6_betti_curve_and_bottleneck() {
    println!("═══════════════════════════════════════════");
    println!("Lesson 6: Betti Curves and Bottleneck Distance");
    println!("═══════════════════════════════════════════");
    println!();

    let mut dg = PersistenceDiagram::new();
    dg.add(0.0, 2.0, 0);
    dg.add(1.0, 3.0, 0);
    dg.add(0.5, 1.5, 1);

    // Betti curve: how many features are alive at each threshold?
    let thresholds = [0.5, 1.0, 1.5, 2.0, 2.5];
    let curve = dg.betti_curve(&thresholds);
    println!("Betti curve (all dimensions):");
    for (t, count) in thresholds.iter().zip(&curve) {
        println!("  t={:.1}: {} alive features", t, count);
    }
    println!();

    // Total persistence.
    let tp1 = dg.total_persistence(1.0);
    let tp2 = dg.total_persistence(2.0);
    println!("Total persistence (p=1): {:.2}", tp1);
    println!("Total persistence (p=2): {:.2}", tp2);
    println!();

    // Bottleneck distance between two diagrams.
    let mut dg2 = PersistenceDiagram::new();
    dg2.add(0.05, 2.05, 0);
    dg2.add(1.05, 3.05, 0);
    dg2.add(0.55, 1.55, 1);
    println!("Bottleneck distance ≈ {:.4}", dg.bottleneck_distance(&dg2));
    println!();
}

// ── Lesson 7: Cellular Sheaves ──────────────────────────────────────────
fn lesson_7_cellular_sheaves() {
    println!("═══════════════════════════════════════════");
    println!("Lesson 7: Cellular Sheaves");
    println!("═══════════════════════════════════════════");
    println!();

    // A constant sheaf: every vertex has stalk R^n, restriction maps are identity.
    let mut cx = SimplicialComplex::new();
    cx.add_edge(0, 1);
    cx.add_edge(1, 2);
    cx.add_edge(2, 3);

    let constant = CellularSheaf::constant(cx.clone(), 3);
    println!("Constant sheaf (stalk_dim=3) on path graph:");
    println!("  stalk dimension: {}", constant.stalk_dimension);
    println!("  global section dim: {}", constant.global_section_dimension());
    println!("  H⁰ dim: {}", constant.cohomology_dimension(0));
    println!("  H¹ dim: {}", constant.cohomology_dimension(1));
    println!();

    // A weighted sheaf: stalks are R, restriction maps scale by weight.
    let weighted = CellularSheaf::from_weights(cx, &[1.0, 0.5, 2.0]);
    println!("Weighted sheaf (weights 1.0, 0.5, 2.0):");
    println!("  stalk dimension: {}", weighted.stalk_dimension);
    println!("  global section dim: {}", weighted.global_section_dimension());
    println!();
}

// ── Lesson 8: Sheaf Laplacian ───────────────────────────────────────────
fn lesson_8_sheaf_laplacian() {
    println!("═══════════════════════════════════════════");
    println!("Lesson 8: Sheaf Laplacian");
    println!("═══════════════════════════════════════════");
    println!();

    // The standard graph Laplacian.
    let edges = vec![(0, 1), (1, 2)];
    let gl = SheafLaplacian::graph_laplacian(3, &edges);
    println!("Graph Laplacian for path 0-1-2 ({}×{}):", gl.dimension, gl.dimension);
    for (i, row) in gl.matrix.iter().enumerate() {
        println!("  row {}: [{:.1}, {:.1}, {:.1}]", i, row[0], row[1], row[2]);
    }
    println!();

    // Multiply by a vector.
    let v = vec![1.0, -1.0, 0.0];
    let lv = gl.mul_vec(&v);
    println!("L × [1, -1, 0] = [{:.1}, {:.1}, {:.1}]", lv[0], lv[1], lv[2]);
    println!();

    // Largest eigenvalue via power iteration.
    let max_eig = gl.largest_eigenvalue(100);
    println!("Largest eigenvalue ≈ {:.4}", max_eig);
    println!();

    // Sheaf Laplacian from a constant sheaf (generalizes graph Laplacian).
    let mut cx = SimplicialComplex::new();
    cx.add_edge(0, 1);
    let sheaf = CellularSheaf::constant(cx, 2);
    let sl = SheafLaplacian::from_sheaf(&sheaf);
    println!("Sheaf Laplacian (2 vertices, stalk_dim=2): {}×{}",
        sl.dimension, sl.dimension);
    println!();
}

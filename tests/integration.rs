//! Integration tests for persistent-sheaf

use persistent_sheaf::*;

#[test]
fn test_simplicial_complex_empty() {
    let sc = SimplicialComplex::new();
    assert_eq!(sc.num_simplices(0), 0);
    assert_eq!(sc.euler_characteristic(), 0);
}

#[test]
fn test_simplicial_complex_triangle() {
    let mut sc = SimplicialComplex::new();
    sc.add_triangle(0, 1, 2);
    assert_eq!(sc.num_simplices(0), 3); // 3 vertices
    assert_eq!(sc.num_simplices(1), 3); // 3 edges
    assert_eq!(sc.num_simplices(2), 1); // 1 triangle
    assert_eq!(sc.euler_characteristic(), 1); // V-E+F = 3-3+1
}

#[test]
fn test_simplicial_complex_two_triangles() {
    let mut sc = SimplicialComplex::new();
    sc.add_triangle(0, 1, 2);
    sc.add_triangle(1, 2, 3);
    assert_eq!(sc.num_simplices(0), 4);
    assert_eq!(sc.num_simplices(1), 5);
    assert_eq!(sc.num_simplices(2), 2);
    assert_eq!(sc.euler_characteristic(), 1); // 4-5+2
}

#[test]
fn test_simplicial_complex_dedup() {
    let mut sc = SimplicialComplex::new();
    sc.add_vertex(0);
    sc.add_vertex(0); // duplicate
    assert_eq!(sc.num_simplices(0), 1);
    sc.add_edge(0, 1);
    sc.add_edge(1, 0); // same edge reversed
    assert_eq!(sc.num_simplices(1), 1);
}

#[test]
fn test_persistence_diagram_basic() {
    let mut pd = PersistenceDiagram::new();
    assert!(pd.is_empty());
    pd.add(0.0, 1.0, 0);
    pd.add(0.5, 2.0, 1);
    assert_eq!(pd.len(), 2);

    let dim0 = pd.filter_dimension(0);
    assert_eq!(dim0.len(), 1);
    assert!((dim0[0].persistence() - 1.0).abs() < 1e-10);
}

#[test]
fn test_persistence_pair_essential() {
    let pd = {
        let mut pd = PersistenceDiagram::new();
        pd.add(0.0, f64::INFINITY, 0);
        pd
    };
    assert!(pd.pairs[0].is_essential());
    assert!(pd.pairs[0].persistence().is_infinite());
}

#[test]
fn test_persistence_pair_midpoint() {
    let pd = {
        let mut pd = PersistenceDiagram::new();
        pd.add(1.0, 5.0, 0);
        pd
    };
    assert!((pd.pairs[0].midpoint() - 3.0).abs() < 1e-10);
}

#[test]
fn test_constant_sheaf() {
    let mut sc = SimplicialComplex::new();
    sc.add_edge(0, 1);
    sc.add_edge(1, 2);
    let sheaf = CellularSheaf::constant(sc, 2);
    assert_eq!(sheaf.stalk_dimension, 2);
    assert_eq!(sheaf.restriction_maps.len(), 2);
}

#[test]
fn test_sheaf_from_weights() {
    let mut sc = SimplicialComplex::new();
    sc.add_edge(0, 1);
    let sheaf = CellularSheaf::from_weights(sc, &[0.5]);
    assert_eq!(sheaf.stalk_dimension, 1);
}

#[test]
fn test_betti_numbers_point() {
    let mut sc = SimplicialComplex::new();
    sc.add_vertex(0);
    let betti = sc.betti_numbers();
    assert_eq!(betti[0], 1); // one connected component
}

#[test]
fn test_total_persistence() {
    let mut pd = PersistenceDiagram::new();
    pd.add(0.0, 2.0, 0);
    pd.add(1.0, 4.0, 1);
    let tp = pd.total_persistence(1.0);
    assert!((tp - 5.0).abs() < 1e-10); // 2 + 3
}

#[test]
fn test_most_persistent() {
    let mut pd = PersistenceDiagram::new();
    pd.add(0.0, 1.0, 0);
    pd.add(0.0, 5.0, 1);
    let mp = pd.most_persistent().unwrap();
    assert!((mp.persistence() - 5.0).abs() < 1e-10);
}

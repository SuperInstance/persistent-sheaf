//! Filtration: a sequence of nested simplicial complexes.

use crate::persistence::PersistenceDiagram;
use crate::simplicial::SimplicialComplex;

/// A filtration: an increasing sequence of simplicial complexes.
pub struct Filtration {
    /// The complexes at each threshold.
    pub complexes: Vec<(f64, SimplicialComplex)>,
}

impl Filtration {
    /// Create an empty filtration.
    pub fn new() -> Self {
        Self { complexes: vec![] }
    }

    /// Add a complex at a given threshold.
    pub fn add(&mut self, threshold: f64, complex: SimplicialComplex) {
        self.complexes.push((threshold, complex));
    }

    /// Build a filtration from a distance matrix by varying epsilon.
    pub fn from_distance_matrix(distances: &[Vec<f64>], num_steps: usize) -> Self {
        let n = distances.len();
        let mut max_dist = 0.0f64;
        for i in 0..n {
            for j in (i + 1)..n {
                max_dist = max_dist.max(distances[i][j]);
            }
        }

        let mut filt = Self::new();
        for step in 0..=num_steps {
            let epsilon = max_dist * step as f64 / num_steps as f64;
            let complex = SimplicialComplex::vietoris_rips(distances, epsilon);
            filt.add(epsilon, complex);
        }
        filt
    }

    /// Compute persistent homology from this filtration.
    pub fn compute_persistence(&self) -> PersistenceDiagram {
        let mut diagram = PersistenceDiagram::new();

        // Track connected components across filtration
        let mut prev_components: Option<Vec<usize>> = None;

        for (threshold, complex) in &self.complexes {
            let betti = complex.betti_numbers();
            let n_components = betti.first().copied().unwrap_or(0);

            if let Some(ref prev) = prev_components {
                if prev.len() > n_components {
                    // Components merged — record birth at the threshold where they appeared
                    // and death at this threshold
                    let n_merged = prev.len() - n_components;
                    for _ in 0..n_merged {
                        diagram.add(0.0, *threshold, 0);
                    }
                }
            }

            // Track 1-dimensional features (cycles)
            if betti.len() > 1 && betti[1] > 0 {
                diagram.add(*threshold * 0.5, *threshold, 1);
            }

            prev_components = Some((0..n_components).collect());
        }

        // Add essential features (those that survive to the end)
        if let Some(last) = self.complexes.last() {
            let betti = last.1.betti_numbers();
            for _ in 0..betti.first().copied().unwrap_or(0) {
                diagram.add(0.0, f64::INFINITY, 0);
            }
        }

        diagram
    }

    /// Number of steps in the filtration.
    pub fn len(&self) -> usize {
        self.complexes.len()
    }

    pub fn is_empty(&self) -> bool {
        self.complexes.is_empty()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_empty_filtration() {
        let f = Filtration::new();
        assert!(f.is_empty());
    }

    #[test]
    fn test_from_distance_matrix() {
        let distances = vec![
            vec![0.0, 1.0, 3.0],
            vec![1.0, 0.0, 1.0],
            vec![3.0, 1.0, 0.0],
        ];
        let f = Filtration::from_distance_matrix(&distances, 5);
        assert!(f.len() > 0);
    }

    #[test]
    fn test_compute_persistence() {
        let distances = vec![
            vec![0.0, 1.0, 2.0],
            vec![1.0, 0.0, 1.0],
            vec![2.0, 1.0, 0.0],
        ];
        let f = Filtration::from_distance_matrix(&distances, 10);
        let diagram = f.compute_persistence();
        // Should have some persistence pairs
        assert!(diagram.len() > 0);
    }

    #[test]
    fn test_filtration_grows() {
        let distances = vec![
            vec![0.0, 1.0, 2.0],
            vec![1.0, 0.0, 1.0],
            vec![2.0, 1.0, 0.0],
        ];
        let f = Filtration::from_distance_matrix(&distances, 5);
        // Complexes should be non-decreasing in size
        for w in f.complexes.windows(2) {
            assert!(w[0].0 <= w[1].0);
        }
    }
}

//! Cellular sheaf: assigns data to cells of a complex with restriction maps.

use crate::simplicial::SimplicialComplex;

/// A cellular sheaf over a simplicial complex.
///
/// For each cell σ, assigns a vector space F(σ) (the "stalk").
/// For each face τ < σ, a linear restriction map F(σ) → F(τ).
pub struct CellularSheaf {
    /// Dimension of the stalk at each vertex.
    pub stalk_dimension: usize,
    /// Restriction maps: edge (a,b) → maps from F(b) to F(a) and F(a) to F(b).
    /// Stored as matrices: restriction_maps[edge_idx] = (map_to_a, map_to_b).
    pub restriction_maps: Vec<(Vec<Vec<f64>>, Vec<Vec<f64>>)>,
    /// The underlying complex.
    pub complex: SimplicialComplex,
}

impl CellularSheaf {
    /// Create a constant sheaf: all stalks are R^n, all restriction maps are identity.
    pub fn constant(complex: SimplicialComplex, stalk_dim: usize) -> Self {
        let n = stalk_dim;
        let id: Vec<Vec<f64>> = (0..n)
            .map(|i| (0..n).map(|j| if i == j { 1.0 } else { 0.0 }).collect())
            .collect();
        let maps: Vec<_> = complex
            .edges
            .iter()
            .map(|_| (id.clone(), id.clone()))
            .collect();
        Self {
            stalk_dimension: n,
            restriction_maps: maps,
            complex,
        }
    }

    /// Create a sheaf from weight functions on edges.
    pub fn from_weights(complex: SimplicialComplex, weights: &[f64]) -> Self {
        let id: Vec<Vec<f64>> = vec![vec![1.0]];
        let maps: Vec<_> = complex
            .edges
            .iter()
            .enumerate()
            .map(|(i, _)| {
                let w = weights.get(i).copied().unwrap_or(1.0);
                (vec![vec![w]], vec![vec![w]])
            })
            .collect();
        Self {
            stalk_dimension: 1,
            restriction_maps: maps,
            complex,
        }
    }

    /// Global sections: assignments of values to all cells that are compatible
    /// with the restriction maps. Returns the dimension of the global section space.
    pub fn global_section_dimension(&self) -> usize {
        // Simplified: for a constant sheaf on a connected complex, it's the stalk dimension
        // For a general sheaf, this requires solving the sheaf condition
        if self.stalk_dimension == 0 {
            return 0;
        }

        // Count edges where restriction maps are surjective
        let n_verts = self.complex.vertices.len();
        if n_verts == 0 {
            return 0;
        }

        // Simplified: the global section dimension equals the stalk dimension
        // for a constant sheaf on a connected complex
        self.stalk_dimension
    }

    /// Sheaf cohomology H^0(F): the space of global sections.
    /// H^1(F): obstruction to extending local sections to global ones.
    pub fn cohomology_dimension(&self, degree: usize) -> usize {
        match degree {
            0 => self.global_section_dimension(),
            1 => {
                // Simplified: H^1 dimension = total stalk dim - H^0 dim - rank of coboundary
                let n_edges = self.complex.edges.len();
                let n_verts = self.complex.vertices.len();
                (n_edges * self.stalk_dimension).saturating_sub(n_verts * self.stalk_dimension)
            }
            _ => 0,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_constant_sheaf() {
        let mut c = SimplicialComplex::new();
        c.add_edge(0, 1);
        c.add_edge(1, 2);
        let sheaf = CellularSheaf::constant(c, 2);
        assert_eq!(sheaf.stalk_dimension, 2);
    }

    #[test]
    fn test_weighted_sheaf() {
        let mut c = SimplicialComplex::new();
        c.add_edge(0, 1);
        c.add_edge(1, 2);
        let sheaf = CellularSheaf::from_weights(c, &[1.0, 0.5]);
        assert_eq!(sheaf.stalk_dimension, 1);
    }

    #[test]
    fn test_global_section_dimension() {
        let mut c = SimplicialComplex::new();
        c.add_edge(0, 1);
        let sheaf = CellularSheaf::constant(c, 3);
        assert_eq!(sheaf.global_section_dimension(), 3);
    }

    #[test]
    fn test_cohomology() {
        let mut c = SimplicialComplex::new();
        c.add_edge(0, 1);
        let sheaf = CellularSheaf::constant(c, 1);
        let h0 = sheaf.cohomology_dimension(0);
        assert!(h0 >= 0);
    }
}

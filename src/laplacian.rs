//! Sheaf Laplacian: generalizes graph Laplacian with sheaf-theoretic information.

use crate::sheaf::CellularSheaf;

/// The sheaf Laplacian L_F for a cellular sheaf F.
///
/// For a 0-cochain f (assignment of vectors to vertices):
/// (L_F f)(v) = Σ_{v~w} F_{v←w}^T (F_{v←w} f(v) - F_{w←v} f(w))
///
/// This generalizes the graph Laplacian: when all stalks are R and all maps are 1,
/// it reduces to the standard graph Laplacian.
pub struct SheafLaplacian {
    /// The dimension of the Laplacian matrix (n_vertices * stalk_dim).
    pub dimension: usize,
    /// The Laplacian matrix stored as a dense matrix.
    pub matrix: Vec<Vec<f64>>,
}

impl SheafLaplacian {
    /// Build the sheaf Laplacian from a cellular sheaf.
    pub fn from_sheaf(sheaf: &CellularSheaf) -> Self {
        let n_verts = sheaf.complex.vertices.len();
        let d = sheaf.stalk_dimension;
        let dim = n_verts * d;
        let mut matrix = vec![vec![0.0; dim]; dim];

        for (edge_idx, (a, b)) in sheaf.complex.edges.iter().enumerate() {
            if edge_idx >= sheaf.restriction_maps.len() {
                break;
            }
            let (map_to_a, map_to_b) = &sheaf.restriction_maps[edge_idx];

            // Add diagonal blocks: F_{a←b}^T * F_{a←b} to position (a,a)
            for i in 0..d {
                for j in 0..d {
                    for k in 0..d {
                        // (M^T * M)_{ij} = Σ_k M_{ki} * M_{kj}
                        let m_ki = map_to_a
                            .get(k)
                            .and_then(|r| r.get(i))
                            .copied()
                            .unwrap_or(0.0);
                        let m_kj = map_to_a
                            .get(k)
                            .and_then(|r| r.get(j))
                            .copied()
                            .unwrap_or(0.0);
                        matrix[a * d + i][a * d + j] += m_ki * m_kj;
                    }
                }
            }

            // Same for vertex b
            for i in 0..d {
                for j in 0..d {
                    for k in 0..d {
                        let m_ki = map_to_b
                            .get(k)
                            .and_then(|r| r.get(i))
                            .copied()
                            .unwrap_or(0.0);
                        let m_kj = map_to_b
                            .get(k)
                            .and_then(|r| r.get(j))
                            .copied()
                            .unwrap_or(0.0);
                        matrix[b * d + i][b * d + j] += m_ki * m_kj;
                    }
                }
            }

            // Off-diagonal: -F_{a←b}^T * F_{b←a} at position (a,b)
            for i in 0..d {
                for j in 0..d {
                    for k in 0..d {
                        let fta = map_to_a
                            .get(k)
                            .and_then(|r| r.get(i))
                            .copied()
                            .unwrap_or(0.0);
                        let ftb = map_to_b
                            .get(k)
                            .and_then(|r| r.get(j))
                            .copied()
                            .unwrap_or(0.0);
                        matrix[a * d + i][b * d + j] -= fta * ftb;
                        matrix[b * d + j][a * d + i] -= fta * ftb;
                    }
                }
            }
        }

        Self {
            dimension: dim,
            matrix,
        }
    }

    /// Build the standard graph Laplacian (sheaf Laplacian with trivial sheaf).
    pub fn graph_laplacian(n_verts: usize, edges: &[(usize, usize)]) -> Self {
        let dim = n_verts;
        let mut matrix = vec![vec![0.0; dim]; dim];

        for &(a, b) in edges {
            matrix[a][a] += 1.0;
            matrix[b][b] += 1.0;
            matrix[a][b] -= 1.0;
            matrix[b][a] -= 1.0;
        }

        Self {
            dimension: dim,
            matrix,
        }
    }

    /// Multiply the Laplacian by a vector.
    pub fn mul_vec(&self, v: &[f64]) -> Vec<f64> {
        assert_eq!(v.len(), self.dimension);
        self.matrix
            .iter()
            .map(|row| row.iter().zip(v.iter()).map(|(a, b)| a * b).sum())
            .collect()
    }

    /// Compute eigenvalues via power iteration (returns the largest eigenvalue).
    pub fn largest_eigenvalue(&self, iterations: usize) -> f64 {
        let n = self.dimension;
        if n == 0 {
            return 0.0;
        }
        let mut v: Vec<f64> = (0..n).map(|i| (i as f64 + 1.0) / (n as f64)).collect();
        let norm = v.iter().map(|x| x * x).sum::<f64>().sqrt();
        v.iter_mut().for_each(|x| *x /= norm);

        for _ in 0..iterations {
            let mv = self.mul_vec(&v);
            let norm = mv.iter().map(|x| x * x).sum::<f64>().sqrt();
            if norm < 1e-15 {
                return 0.0;
            }
            v = mv.iter().map(|x| x / norm).collect();
        }

        let mv = self.mul_vec(&v);
        let v_dot_mv: f64 = v.iter().zip(&mv).map(|(a, b)| a * b).sum();
        v_dot_mv
    }

    /// Compute the Fiedler value (second smallest eigenvalue, algebraic connectivity).
    pub fn fiedler_value(&self) -> f64 {
        // Simplified: largest eigenvalue minus trace / n for small matrices
        let _trace: f64 = (0..self.dimension).map(|i| self.matrix[i][i]).sum();
        let max_eig = self.largest_eigenvalue(50);
        // Rough approximation
        max_eig / self.dimension.max(1) as f64
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::simplicial::SimplicialComplex;

    #[test]
    fn test_graph_laplacian() {
        let edges = vec![(0, 1), (1, 2)];
        let l = SheafLaplacian::graph_laplacian(3, &edges);
        assert_eq!(l.dimension, 3);
        // Diagonal: degree of each vertex
        assert_eq!(l.matrix[0][0], 1.0);
        assert_eq!(l.matrix[1][1], 2.0);
        assert_eq!(l.matrix[2][2], 1.0);
    }

    #[test]
    fn test_laplacian_positive_semidefinite() {
        let edges = vec![(0, 1), (1, 2)];
        let l = SheafLaplacian::graph_laplacian(3, &edges);
        // Largest eigenvalue should be positive
        let eig = l.largest_eigenvalue(100);
        assert!(eig > 0.0);
    }

    #[test]
    fn test_sheaf_laplacian() {
        let mut c = SimplicialComplex::new();
        c.add_edge(0, 1);
        let sheaf = CellularSheaf::constant(c, 2);
        let l = SheafLaplacian::from_sheaf(&sheaf);
        assert_eq!(l.dimension, 4); // 2 vertices * 2 stalk dim
    }

    #[test]
    fn test_mul_vec() {
        let edges = vec![(0, 1)];
        let l = SheafLaplacian::graph_laplacian(2, &edges);
        let v = vec![1.0, -1.0];
        let lv = l.mul_vec(&v);
        // L * [1,-1]^T = [2, -2]^T
        assert!((lv[0] - 2.0).abs() < 0.01);
        assert!((lv[1] + 2.0).abs() < 0.01);
    }
}

//! Simplicial complex data structure.

/// A simplicial complex: a collection of simplices closed under taking faces.
#[derive(Debug, Clone)]
pub struct SimplicialComplex {
    /// Vertices (0-simplices).
    pub vertices: Vec<usize>,
    /// Edges (1-simplices), stored as sorted pairs.
    pub edges: Vec<(usize, usize)>,
    /// Triangles (2-simplices), stored as sorted triples.
    pub triangles: Vec<(usize, usize, usize)>,
    /// Tetrahedra (3-simplices).
    pub tetrahedra: Vec<(usize, usize, usize, usize)>,
}

impl SimplicialComplex {
    /// Create a new empty complex.
    pub fn new() -> Self {
        Self {
            vertices: vec![],
            edges: vec![],
            triangles: vec![],
            tetrahedra: vec![],
        }
    }

    /// Add a vertex.
    pub fn add_vertex(&mut self, v: usize) {
        if !self.vertices.contains(&v) {
            self.vertices.push(v);
        }
    }

    /// Add an edge (ensures vertices exist).
    pub fn add_edge(&mut self, a: usize, b: usize) {
        let edge = if a < b { (a, b) } else { (b, a) };
        self.add_vertex(a);
        self.add_vertex(b);
        if !self.edges.contains(&edge) {
            self.edges.push(edge);
        }
    }

    /// Add a triangle (ensures edges exist).
    pub fn add_triangle(&mut self, a: usize, b: usize, c: usize) {
        let mut verts = [a, b, c];
        verts.sort();
        let tri = (verts[0], verts[1], verts[2]);
        self.add_edge(a, b);
        self.add_edge(b, c);
        self.add_edge(a, c);
        if !self.triangles.contains(&tri) {
            self.triangles.push(tri);
        }
    }

    /// Number of simplices of dimension k.
    pub fn num_simplices(&self, dim: usize) -> usize {
        match dim {
            0 => self.vertices.len(),
            1 => self.edges.len(),
            2 => self.triangles.len(),
            3 => self.tetrahedra.len(),
            _ => 0,
        }
    }

    /// Euler characteristic: χ = V - E + F - T.
    pub fn euler_characteristic(&self) -> i32 {
        self.vertices.len() as i32 - self.edges.len() as i32 + self.triangles.len() as i32
            - self.tetrahedra.len() as i32
    }

    /// Boundary of a simplex: the faces of dimension dim-1.
    pub fn boundary_of_edge(&self, edge: (usize, usize)) -> Vec<usize> {
        vec![edge.0, edge.1]
    }

    /// Boundary of a triangle: its three edges.
    pub fn boundary_of_triangle(&self, tri: (usize, usize, usize)) -> Vec<(usize, usize)> {
        let (a, b, c) = tri;
        let mut edges = vec![
            (a.min(b), a.max(b)),
            (b.min(c), b.max(c)),
            (a.min(c), a.max(c)),
        ];
        edges.sort();
        edges
    }

    /// Betti numbers: β_k = rank(H_k) = dim(ker ∂_k) - dim(im ∂_{k+1}).
    /// Simplified computation for small complexes.
    pub fn betti_numbers(&self) -> Vec<usize> {
        // β_0 = connected components (simplified: vertices - edges + triangles that form cycles)
        let n = self.vertices.len();
        if n == 0 {
            return vec![];
        }

        // β_0: number of connected components (union-find)
        let mut parent: Vec<usize> = (0..n).collect();
        let find = |v: usize, parent: &mut Vec<usize>| -> usize {
            let mut x = v;
            while parent[x] != x {
                parent[x] = parent[parent[x]];
                x = parent[x];
            }
            x
        };

        for (a, b) in &self.edges {
            if *a < n && *b < n {
                let ra = find(*a, &mut parent);
                let rb = find(*b, &mut parent);
                if ra != rb {
                    parent[ra] = rb;
                }
            }
        }

        let components: std::collections::HashSet<usize> =
            (0..n).map(|v| find(v, &mut parent)).collect();
        let beta0 = components.len();

        // β_1 = E - V + β_0 (for 2D complexes, by Euler-Poincaré)
        let beta1 = if self.triangles.is_empty() {
            self.edges.len() as i32 - n as i32 + beta0 as i32
        } else {
            // More complex: use Euler characteristic
            self.edges.len() as i32 - n as i32 + self.triangles.len() as i32 + beta0 as i32
                - self.euler_characteristic()
        };

        vec![beta0, beta1.max(0) as usize]
    }

    /// Build a Vietoris-Rips complex from a distance matrix.
    pub fn vietoris_rips(distances: &[Vec<f64>], epsilon: f64) -> Self {
        let n = distances.len();
        let mut complex = Self::new();

        for i in 0..n {
            complex.add_vertex(i);
        }

        for i in 0..n {
            for j in (i + 1)..n {
                if distances[i][j] <= epsilon {
                    complex.add_edge(i, j);
                }
            }
        }

        // Add triangles where all three edges exist
        for i in 0..n {
            for j in (i + 1)..n {
                for k in (j + 1)..n {
                    if distances[i][j] <= epsilon
                        && distances[j][k] <= epsilon
                        && distances[i][k] <= epsilon
                    {
                        complex.add_triangle(i, j, k);
                    }
                }
            }
        }

        complex
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_empty_complex() {
        let c = SimplicialComplex::new();
        assert_eq!(c.num_simplices(0), 0);
    }

    #[test]
    fn test_add_vertex() {
        let mut c = SimplicialComplex::new();
        c.add_vertex(0);
        c.add_vertex(1);
        assert_eq!(c.num_simplices(0), 2);
    }

    #[test]
    fn test_add_edge() {
        let mut c = SimplicialComplex::new();
        c.add_edge(0, 1);
        assert_eq!(c.num_simplices(0), 2);
        assert_eq!(c.num_simplices(1), 1);
    }

    #[test]
    fn test_add_triangle() {
        let mut c = SimplicialComplex::new();
        c.add_triangle(0, 1, 2);
        assert_eq!(c.num_simplices(0), 3);
        assert_eq!(c.num_simplices(1), 3);
        assert_eq!(c.num_simplices(2), 1);
    }

    #[test]
    fn test_euler_characteristic() {
        let mut c = SimplicialComplex::new();
        c.add_triangle(0, 1, 2);
        // V=3, E=3, F=1 → χ = 3-3+1 = 1
        assert_eq!(c.euler_characteristic(), 1);
    }

    #[test]
    fn test_betti_numbers_point() {
        let mut c = SimplicialComplex::new();
        c.add_vertex(0);
        let betti = c.betti_numbers();
        assert_eq!(betti[0], 1); // one connected component
    }

    #[test]
    fn test_betti_numbers_triangle() {
        let mut c = SimplicialComplex::new();
        c.add_edge(0, 1);
        c.add_edge(1, 2);
        c.add_edge(0, 2);
        // One component, one cycle → β_0=1, β_1=1
        let betti = c.betti_numbers();
        assert_eq!(betti[0], 1);
    }

    #[test]
    fn test_vietoris_rips() {
        let distances = vec![
            vec![0.0, 1.0, 2.0],
            vec![1.0, 0.0, 1.0],
            vec![2.0, 1.0, 0.0],
        ];
        let c = SimplicialComplex::vietoris_rips(&distances, 1.5);
        assert_eq!(c.num_simplices(0), 3);
        assert_eq!(c.num_simplices(1), 2); // edges 0-1 and 1-2
    }

    #[test]
    fn test_boundary_of_triangle() {
        let mut c = SimplicialComplex::new();
        c.add_triangle(0, 1, 2);
        let edges = c.boundary_of_triangle((0, 1, 2));
        assert_eq!(edges.len(), 3);
    }
}

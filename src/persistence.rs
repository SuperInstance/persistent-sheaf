//! Persistence diagrams for tracking topological features across scales.

/// A birth-death pair in a persistence diagram.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PersistencePair {
    pub birth: f64,
    pub death: f64,
    pub dimension: usize,
}

impl PersistencePair {
    /// Create a new persistence pair.
    pub fn new(birth: f64, death: f64, dimension: usize) -> Self {
        Self {
            birth,
            death,
            dimension,
        }
    }

    /// Persistence: how long the feature survives.
    pub fn persistence(&self) -> f64 {
        self.death - self.birth
    }

    /// Whether this is an essential feature (never dies).
    pub fn is_essential(&self) -> bool {
        self.death == f64::INFINITY
    }

    /// Midpoint of the birth-death interval.
    pub fn midpoint(&self) -> f64 {
        (self.birth + self.death) / 2.0
    }
}

/// A persistence diagram: collection of birth-death pairs.
#[derive(Debug, Clone)]
pub struct PersistenceDiagram {
    pub pairs: Vec<PersistencePair>,
}

impl PersistenceDiagram {
    /// Create an empty diagram.
    pub fn new() -> Self {
        Self { pairs: vec![] }
    }

    /// Add a pair.
    pub fn add(&mut self, birth: f64, death: f64, dimension: usize) {
        self.pairs
            .push(PersistencePair::new(birth, death, dimension));
    }

    /// Filter pairs by dimension.
    pub fn filter_dimension(&self, dim: usize) -> Vec<&PersistencePair> {
        self.pairs.iter().filter(|p| p.dimension == dim).collect()
    }

    /// Number of pairs.
    pub fn len(&self) -> usize {
        self.pairs.len()
    }

    pub fn is_empty(&self) -> bool {
        self.pairs.is_empty()
    }

    /// Bottleneck distance between two diagrams.
    pub fn bottleneck_distance(&self, other: &Self) -> f64 {
        // Simplified: max of minimum distances between all pairs
        let mut max_min = 0.0f64;
        for p1 in &self.pairs {
            let min_dist = other
                .pairs
                .iter()
                .map(|p2| ((p1.birth - p2.birth).abs()).max((p1.death - p2.death).abs()))
                .fold(f64::INFINITY, f64::min);
            max_min = max_min.max(min_dist);
        }
        max_min
    }

    /// Total persistence: sum of all persistence values.
    pub fn total_persistence(&self, power: f64) -> f64 {
        self.pairs.iter().map(|p| p.persistence().powf(power)).sum()
    }

    /// The most persistent feature.
    pub fn most_persistent(&self) -> Option<&PersistencePair> {
        self.pairs
            .iter()
            .max_by(|a, b| a.persistence().partial_cmp(&b.persistence()).unwrap())
    }

    /// Betti curve: number of alive features at each threshold.
    pub fn betti_curve(&self, thresholds: &[f64]) -> Vec<usize> {
        thresholds
            .iter()
            .map(|&t| {
                self.pairs
                    .iter()
                    .filter(|p| p.birth <= t && p.death > t)
                    .count()
            })
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_persistence_pair() {
        let p = PersistencePair::new(0.0, 1.0, 0);
        assert_eq!(p.persistence(), 1.0);
        assert!(!p.is_essential());
        assert_eq!(p.midpoint(), 0.5);
    }

    #[test]
    fn test_essential_pair() {
        let p = PersistencePair::new(0.0, f64::INFINITY, 0);
        assert!(p.is_essential());
        assert_eq!(p.persistence(), f64::INFINITY);
    }

    #[test]
    fn test_diagram_filter() {
        let mut d = PersistenceDiagram::new();
        d.add(0.0, 1.0, 0);
        d.add(0.5, 2.0, 1);
        d.add(0.3, 0.8, 0);
        assert_eq!(d.filter_dimension(0).len(), 2);
        assert_eq!(d.filter_dimension(1).len(), 1);
    }

    #[test]
    fn test_most_persistent() {
        let mut d = PersistenceDiagram::new();
        d.add(0.0, 1.0, 0);
        d.add(0.0, 3.0, 0);
        d.add(0.5, 0.7, 0);
        let mp = d.most_persistent().unwrap();
        assert_eq!(mp.persistence(), 3.0);
    }

    #[test]
    fn test_total_persistence() {
        let mut d = PersistenceDiagram::new();
        d.add(0.0, 1.0, 0);
        d.add(0.0, 2.0, 0);
        assert_eq!(d.total_persistence(1.0), 3.0);
    }

    #[test]
    fn test_betti_curve() {
        let mut d = PersistenceDiagram::new();
        d.add(0.0, 2.0, 0);
        d.add(1.0, 3.0, 0);
        let curve = d.betti_curve(&[0.5, 1.5, 2.5]);
        assert_eq!(curve[0], 1); // only first feature alive at 0.5
        assert_eq!(curve[1], 2); // both alive at 1.5
        assert_eq!(curve[2], 1); // only second alive at 2.5
    }

    #[test]
    fn test_bottleneck_distance() {
        let mut d1 = PersistenceDiagram::new();
        d1.add(0.0, 1.0, 0);
        let mut d2 = PersistenceDiagram::new();
        d2.add(0.1, 1.1, 0);
        let dist = d1.bottleneck_distance(&d2);
        assert!(dist < 0.2);
    }
}

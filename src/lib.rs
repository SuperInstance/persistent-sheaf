//! # Persistent Sheaf
//!
//! Topological data analysis via persistent sheaf cohomology.
//! Combines persistent homology with sheaf theory for multi-modal data fusion.
//!
//! # Key Concepts
//!
//! - **Cellular sheaf**: Assigns data to cells of a simplicial complex with
//!   restriction maps between them
//! - **Sheaf Laplacian**: Generalizes the graph Laplacian to encode both
//!   geometric and non-geometric information
//! - **Persistence**: Tracks how topological features appear and disappear
//!   across scale parameters

mod filtration;
mod laplacian;
mod persistence;
mod sheaf;
mod simplicial;

pub use filtration::Filtration;
pub use laplacian::SheafLaplacian;
pub use persistence::PersistenceDiagram;
pub use sheaf::CellularSheaf;
pub use simplicial::SimplicialComplex;

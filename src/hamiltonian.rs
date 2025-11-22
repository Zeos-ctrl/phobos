use crate::Complex;
use crate::Matrix2;

use std::f64::consts::PI;

/// Simulates Hamiltonian time evolution e^(-iHt) for a single qubit.
/// 
/// For a Hermitian matrix H (the Hamiltonian), this computes the unitary operator
/// U(t) = e^(-iHt) that represents quantum time evolution for duration t.
/// 
/// # Mathematical Background
/// 
/// If H has eigendecomposition H = Σᵢ λᵢ|vᵢ⟩⟨vᵢ|, then:
/// e^(-iHt) = Σᵢ e^(-iλᵢt)|vᵢ⟩⟨vᵢ|
/// 
/// Each eigenvector |vᵢ⟩ picks up a phase e^(-iλᵢt). We store these as
/// (θ, |v⟩⟨v|) pairs where θ = -λt/π, so the phase becomes e^(iπθ).
/// This representation works naturally with quantum phase estimation.
/// 
/// # Example
/// 
/// ```rust
/// use nalgebra::{Matrix2, Complex};
/// use phobos::HamiltonianSimulation;
/// 
/// // Pauli-Z Hamiltonian: H = [[1, 0], [0, -1]]
/// let h = Matrix2::new(
///     Complex::new(1.0, 0.0), Complex::new(0.0, 0.0),
///     Complex::new(0.0, 0.0), Complex::new(-1.0, 0.0)
/// );
/// 
/// let ham_sim = HamiltonianSimulation::new(h, 0.5);
/// let unitary = ham_sim.get_matrix(); // Returns e^(-iHt)
/// ```
pub struct HamiltonianSimulation {
    /// The Hamiltonian matrix H (must be Hermitian)
    h: Matrix2,
    
    /// Time duration for evolution
    t: f64,
    
    /// Exponent for fractional powers of the gate (1.0 = full gate)
    /// Allows computing U^α for any α
    exponent: f64,
    
    /// Eigendecomposition stored as (θ, projector) pairs
    /// where θ = -λt/π and projector = |v⟩⟨v|
    /// This allows reconstruction: U = Σᵢ e^(iπθᵢ)|vᵢ⟩⟨vᵢ|
    eigen_components: Vec<(f64, Matrix2)>,
}

impl HamiltonianSimulation {
    /// Creates a new Hamiltonian simulation with exponent = 1.0
    /// 
    /// This is a convenience constructor that creates the full unitary U = e^(-iHt).
    /// For fractional powers (like U^(1/2) or U^2), use `with_exponent()` instead.
    /// 
    /// # Arguments
    /// 
    /// * `h` - The Hamiltonian matrix (must be Hermitian: H = H†)
    /// * `t` - Time duration for evolution
    /// 
    /// # Example
    /// 
    /// ```rust
    /// let h = Matrix2::new(/* Pauli-Z */);
    /// let ham_sim = HamiltonianSimulation::new(h, 0.5);
    /// ```
    pub fn new(h: Matrix2, t: f64) -> Self {
        Self::with_exponent(h, t, 1.0)
    }
    
    /// Creates a Hamiltonian simulation with a custom exponent
    /// 
    /// Computes U^α = e^(-iHtα) by eigendecomposing H and storing the components
    /// in a form suitable for quantum phase estimation.
    /// 
    /// # Arguments
    /// 
    /// * `h` - The Hamiltonian matrix (must be Hermitian)
    /// * `t` - Time duration for evolution
    /// * `exponent` - Power to raise the gate to (α). Use 1.0 for standard evolution,
    ///                0.5 for square root, 2.0 for squared gate, etc.
    /// 
    /// # Algorithm
    /// 
    /// 1. Eigendecompose H = Σᵢ λᵢ|vᵢ⟩⟨vᵢ|
    /// 2. For each eigenvalue λᵢ, compute θᵢ = -λᵢt/π
    /// 3. Store (θᵢ, |vᵢ⟩⟨vᵢ|) pairs as eigen components
    /// 4. Later, reconstruct U^α = Σᵢ e^(iπθᵢα)|vᵢ⟩⟨vᵢ|
    /// 
    /// # Why θ = -λt/π?
    /// 
    /// We want e^(-iλt), which we can write as e^(iπθ) where θ = -λt/π.
    /// This θ representation is convenient for quantum phase estimation, which
    /// estimates phases in units of π.
    pub fn with_exponent(h: Matrix2, t: f64, exponent: f64) -> Self {
        let (lambda1, lambda2) = h.eigenvalues_hermitian();
        let (v1, v2) = h.eigenvectors_hermitian(lambda1, lambda2);
        
        let mut eigen_components = Vec::new();
        
        // For first eigenvector/eigenvalue pair
        let theta1 = -lambda1 * t / PI;
        let proj1 = v1.outer_product();
        eigen_components.push((theta1, proj1));
        
        let theta2 = -lambda2 * t / PI;
        let proj2 = v2.outer_product();
        eigen_components.push((theta2, proj2));
        
        Self { h, t, exponent, eigen_components }
    }
    
    /// Reconstructs the full unitary matrix U^α = e^(-iHtα)
    /// 
    /// Uses the spectral decomposition to rebuild the matrix:
    /// U^α = Σᵢ e^(iπθᵢα)|vᵢ⟩⟨vᵢ|
    /// 
    /// Each eigenvector contributes its projector weighted by its phase.
    /// 
    /// # Returns
    /// 
    /// The 2×2 unitary matrix representing the time evolution operator
    /// 
    /// # Example
    /// 
    /// ```rust
    /// let ham_sim = HamiltonianSimulation::new(h, t);
    /// let u_matrix = ham_sim.get_matrix();
    /// 
    /// // Use in a quantum circuit
    /// circuit.add_gate(Gate::U { target: 0, matrix: u_matrix });
    /// ```
    pub fn get_matrix(&self) -> Matrix2 {
        // Reconstruct U = sum_i e^(i*pi*theta*exponent) * |v_i><v_i|
        let mut result = Matrix2::zeros();
        
        for (theta, projector) in &self.eigen_components {
            // Compute phase: e^(iπθα)
            let phase = Complex::new(0.0, PI * theta * self.exponent).exp();
            
            // Add this eigenspace's contribution
            let scaled_projection = projector.scale(&phase);
            result = result.add(&scaled_projection);
        }
        
        result
    }
    
    /// Returns a reference to the stored eigen components
    /// 
    /// Each component is a tuple (θ, P) where:
    /// - θ is the phase angle (in units of π)
    /// - P is the projector |v⟩⟨v| onto that eigenvector
    /// 
    /// This is useful for algorithms that need direct access to the spectral
    /// decomposition, such as when implementing controlled versions of this gate.
    pub fn eigen_components(&self) -> &Vec<(f64, Matrix2)> {
        &self.eigen_components
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_hamiltonian_simulation() {
        let h = Matrix2::new(
            Complex::new(1.0, 0.0), Complex::new(0.0, 0.0),
            Complex::new(0.0, 0.0), Complex::new(-1.0, 0.0)
        );
        
        let ham_sim = HamiltonianSimulation::new(h, PI/4.0);
        let u = ham_sim.get_matrix();
        
        // Test: U × U† should equal identity matrix
        let u_dagger = u.conjugate_transpose();
        let product = u.multiply(&u_dagger);
        
        // Check that product is identity
        assert!((product[(0,0)].real - 1.0).abs() < 1e-10);
        assert!(product[(0,0)].imag.abs() < 1e-10);
        assert!(product[(0,1)].real.abs() < 1e-10);
        assert!(product[(0,1)].imag.abs() < 1e-10);
        assert!(product[(1,0)].real.abs() < 1e-10);
        assert!(product[(1,0)].imag.abs() < 1e-10);
        assert!((product[(1,1)].real - 1.0).abs() < 1e-10);
        assert!(product[(1,1)].imag.abs() < 1e-10);
    }
}

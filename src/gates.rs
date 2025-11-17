use crate::state::QuantumState;
use crate::Complex;

use rand::Rng;

/// Applies a Hadamard gate to the specified qubit.
///
/// The Hadamard gate creates superposition, transforming:
/// - |0⟩ → (1/√2)(|0⟩ + |1⟩)
/// - |1⟩ → (1/√2)(|0⟩ - |1⟩)
///
/// # Arguments
/// * `state` - The quantum state to modify
/// * `target_qubit` - Index of the qubit to apply the gate to (0-indexed)
///
/// # Examples
/// ```
/// use phobos::{QuantumState, gates};
/// let mut state = QuantumState::new(1);
/// gates::apply_hadamard(&mut state, 0);
/// ```
pub fn apply_hadamard(state: &mut QuantumState, target_qubit: usize) {
    let num_amplitudes = state.amplitudes.len();
    
    for i in 0..num_amplitudes {
        // Check if bit at position target_qubit in i is 0 using a bitwise mask
        if (i & (1 << target_qubit)) == 0 {
            let j = i + (1 << target_qubit);  // Same as i + 2^target_qubit

            let old_i = state.amplitudes[i];
            let old_j = state.amplitudes[j];

            let factor = 1.0 / f64::sqrt(2.0);

            // Compute new values
            let new_i = old_i.add(&old_j).scale(factor);
            let new_j = old_i.subtract(&old_j).scale(factor);

            // Update the state
            state.amplitudes[i] = new_i;
            state.amplitudes[j] = new_j;
            
        }
    } 
}

/// Applies a CNOT gate to the specified qubit.
///
/// The CNOT gate flips the target qubit if the control bit is |1⟩:
///
/// If our control bit is 0 and the target bit is 1 we have this truth
/// table;
///
/// - |00> -> |00>
/// - |01> -> |01>
/// - |10> -> |11>
/// - |11> -> |10>
///
/// # Arguments
/// * `state` - The quantum state to modify
/// * `control_qubit` - Index of the qubit to apply the gate to (0-indexed)
/// * `target_qubit` - Index of the qubit to apply the gate to (0-indexed)
///
/// # Examples
/// ```
/// use phobos::{QuantumState, gates};
/// let mut state = QuantumState::new(2);
/// gates::apply_cnot(&mut state, 0, 1);
/// ```
pub fn apply_cnot(state: &mut QuantumState, control_qubit: usize, target_qubit: usize) {
    let num_amplitudes = state.amplitudes.len();

    for i in 0..num_amplitudes {
        // First check: is the control bit set to 1?
        if (i & (1 << control_qubit)) != 0 {
            // only process when target bit is 0
            if (i & (1 << target_qubit)) == 0 {
                // Calculate the swap partner
                let j = i ^ (1 << target_qubit);  // XOR flips the target bit
                
                // Swap amplitudes[i] and amplitudes[j]
                let temp = state.amplitudes[i];
                state.amplitudes[i] = state.amplitudes[j];
                state.amplitudes[j] = temp;
            }
        }
    }
}

/// Applies a Controlled-Z gate to the specified qubits.
///
/// The Controlled-Z gate applies a phase flip to the target qubit if and 
/// only if the control qubit is in state |1⟩. This results in adding a 
/// negative phase only to the |11⟩ basis state:
///
/// - |00⟩ → |00⟩
/// - |01⟩ → |01⟩
/// - |10⟩ → |10⟩
/// - |11⟩ → -|11⟩
///
/// # Arguments
/// * `state` - The quantum state to modify
/// * `control_qubit` - Index of the control qubit (0-indexed)
/// * `target_qubit` - Index of the target qubit (0-indexed)
///
/// # Examples
/// ```
/// use phobos::{QuantumState, gates};
/// let mut state = QuantumState::new(2);
/// gates::apply_cz(&mut state, 0, 1);
/// ```
pub fn apply_cz(state: &mut QuantumState, control_qubit: usize, target_qubit: usize) {
    let num_amplitudes = state.amplitudes.len();
    
    for i in 0..num_amplitudes {
        // Check if BOTH control and target bits are 1
        let control_mask = 1 << control_qubit;
        let target_mask = 1 << target_qubit;
        
        if (i & control_mask) != 0 && (i & target_mask) != 0 {
            // Both qubits are |1⟩, so add negative phase
            state.amplitudes[i] = state.amplitudes[i].scale(-1.0);
        }
    }
}

/// Applies a Pauli-Z gate to the specified qubit.
///
/// The Z gate applies a phase flip, leaving |0⟩ unchanged but adding
/// a negative phase to |1⟩:
/// - |0⟩ → |0⟩
/// - |1⟩ → -|1⟩
///
/// # Arguments
/// * `state` - The quantum state to modify
/// * `target_qubit` - Index of the qubit to apply the gate to (0-indexed)
///
/// # Examples
/// ```
/// use phobos::{QuantumState, gates};
/// let mut state = QuantumState::new(1);
/// gates::apply_z(&mut state, 0);
/// ```
pub fn apply_z(state: &mut QuantumState, target_qubit: usize) {
    let num_amplitudes = state.amplitudes.len();
    
    for i in 0..num_amplitudes {
        // Check if bit at position target_qubit in i is 0 using a bitwise mask
        if (i & (1 << target_qubit)) == 0 {
            let j = i + (1 << target_qubit);  // Same as i + 2^target_qubit

            let old_j = state.amplitudes[j];
            let new_j = old_j.scale(-1.0);

            // Update the state
            state.amplitudes[j] = new_j;
        }
    }
}

/// Applies a Pauli-X gate to the specified qubit.
///
/// The X gate is the quantum NOT gate, flipping the qubit state:
/// - |0⟩ → |1⟩
/// - |1⟩ → |0⟩
///
/// # Arguments
/// * `state` - The quantum state to modify
/// * `target_qubit` - Index of the qubit to apply the gate to (0-indexed)
///
/// # Examples
/// ```
/// use phobos::{QuantumState, gates};
/// let mut state = QuantumState::new(1);
/// gates::apply_x(&mut state, 0);
/// ```
pub fn apply_x(state: &mut QuantumState, target_qubit: usize) {
    let num_amplitudes = state.amplitudes.len();
    
    for i in 0..num_amplitudes {
        // Check if bit at position target_qubit in i is 0 using a bitwise mask
        if (i & (1 << target_qubit)) == 0 {
            let j = i + (1 << target_qubit);  // Same as i + 2^target_qubit

            let old_i = state.amplitudes[i]; // Alpha
            let old_j = state.amplitudes[j]; // Beta

            // Swap alpha and beta
            let new_i = old_j;
            let new_j = old_i;

            // Update the state
            state.amplitudes[i] = new_i;
            state.amplitudes[j] = new_j;
        }
    }
}

/// Applies a Pauli-Y gate to the specified qubit.
///
/// The Y gate combines a bit flip with a phase change:
/// - |0⟩ → i|1⟩
/// - |1⟩ → -i|0⟩
///
/// # Arguments
/// * `state` - The quantum state to modify
/// * `target_qubit` - Index of the qubit to apply the gate to (0-indexed)
///
/// # Examples
/// ```
/// use phobos::{QuantumState, gates};
/// let mut state = QuantumState::new(1);
/// gates::apply_y(&mut state, 0);
/// ```
pub fn apply_y(state: &mut QuantumState, target_qubit: usize) {
    let num_amplitudes = state.amplitudes.len();
    let imag = Complex::new(0.0, 1.0);
    
    for i in 0..num_amplitudes {
        // Check if bit at position target_qubit in i is 0 using a bitwise mask
        if (i & (1 << target_qubit)) == 0 {
            let j = i + (1 << target_qubit);  // Same as i + 2^target_qubit

            let old_i = state.amplitudes[i]; // Alpha
            let old_j = state.amplitudes[j]; // Beta

            // Swap alpha and beta
            let new_i = old_j.multiply(&imag.conjugate());
            let new_j = old_i.multiply(&imag);

            // Update the state
            state.amplitudes[i] = new_i;
            state.amplitudes[j] = new_j;
        }
    }
}

/// Applies an Identity gate to the specified qubit.
///
/// The Identity gate performs no operation, leaving the state unchanged.
/// This is primarily used as a placeholder in quantum circuits.
///
/// # Arguments
/// * `state` - The quantum state (unchanged)
/// * `target_qubit` - Index of the qubit (unused)
pub fn apply_identity(state: &mut QuantumState, target_qubit: usize) {
    let num_amplitudes = state.amplitudes.len();

    for i in 0..num_amplitudes {
        if (i & (1 << target_qubit)) == 0 {
            let j = i + (1 << target_qubit);

            let identity_i = state.amplitudes[i];
            let identity_j = state.amplitudes[j];

            state.amplitudes[i] = identity_i;
            state.amplitudes[j] = identity_j;
        }
    }
}

/// Measures a specific qubit and collapses the quantum state.
///
/// Unlike [`QuantumState::measure`] which measures all qubits, this function
/// measures only the specified qubit, leaving other qubits in superposition.
/// The state is collapsed based on the measurement outcome and renormalized.
///
/// # Arguments
/// * `state` - The quantum state to measure and collapse
/// * `target_qubit` - Index of the qubit to measure (0-indexed)
///
/// # Returns
/// The measurement result: 0 or 1
///
/// # Examples
/// ```
/// use phobos::{QuantumState, gates};
/// 
/// // Create superposition on qubit 0
/// let mut state = QuantumState::new(2);
/// gates::apply_hadamard(&mut state, 0);
/// 
/// // Measure only qubit 0 (qubit 1 remains in |0⟩)
/// let result = gates::measure_qubit(&mut state, 0);
/// // result is either 0 or 1, and the state has collapsed accordingly
/// ```
pub fn measure_qubit(state: &mut QuantumState, target_qubit: usize) -> u8 {
    let mut prob_zero = 0.0;
    
    // Get the probabiities where the target qubit is 0
    for i in 0..state.amplitudes.len() {
        if (i & (1 << target_qubit)) == 0 {
            prob_zero += state.amplitudes[i].magnitude_squared();
        }
    }
    
    // Measure one probability
    let mut rng = rand::rng();
    let r: f64 = rng.random();

    let outcome = if r < prob_zero { 0 } else { 1 };

    for i in 0..state.amplitudes.len() {
        let bit_value = if (i & (1 << target_qubit)) == 0 { 0 } else { 1 };
        if bit_value != outcome {
            state.amplitudes[i] = Complex::new(0.0, 0.0);
        }
    }
    state.normalize();
    return outcome;
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::state::QuantumState;

    #[test]
    fn test_hadamard_on_zero_state() {
        let mut state = QuantumState::new(1);
        apply_hadamard(&mut state, 0);
        
        // After H|0⟩, should get (1/√2)|0⟩ + (1/√2)|1⟩
        let expected = 1.0 / f64::sqrt(2.0);
        assert!((state.amplitudes[0].real - expected).abs() < 1e-10);
        assert!((state.amplitudes[1].real - expected).abs() < 1e-10);
        assert!(state.amplitudes[0].imag.abs() < 1e-10);
        assert!(state.amplitudes[1].imag.abs() < 1e-10);
    }

    #[test]
    fn test_cnot_creates_bell_state() {
        let mut state = QuantumState::new(2);
        
        // Apply H to qubit 0: |00⟩ → (1/√2)(|00⟩ + |10⟩)
        apply_hadamard(&mut state, 0);
        
        // Apply CNOT: (1/√2)(|00⟩ + |10⟩) → (1/√2)(|00⟩ + |11⟩)
        apply_cnot(&mut state, 0, 1);
        
        let expected = 1.0 / f64::sqrt(2.0);
        
        // Should have amplitude 1/√2 for |00⟩ (index 0)
        assert!((state.amplitudes[0].real - expected).abs() < 1e-10);
        
        // Should have amplitude 0 for |01⟩ (index 1) and |10⟩ (index 2)
        assert!(state.amplitudes[1].magnitude_squared() < 1e-10);
        assert!(state.amplitudes[2].magnitude_squared() < 1e-10);
        
        // Should have amplitude 1/√2 for |11⟩ (index 3)
        assert!((state.amplitudes[3].real - expected).abs() < 1e-10);
    }

    #[test]
    fn test_z_gate() {
        use crate::state::QuantumState;
        use crate::complex::Complex;
        
        // Create a state in superposition: (1/√2)|0⟩ + (1/√2)|1⟩
        let mut state = QuantumState::new(1);
        state.amplitudes[0] = Complex::new(1.0 / f64::sqrt(2.0), 0.0);
        state.amplitudes[1] = Complex::new(1.0 / f64::sqrt(2.0), 0.0);
        
        // Apply Z gate
        apply_z(&mut state, 0);
        
        // After Z: (1/√2)|0⟩ - (1/√2)|1⟩
        // First amplitude should stay positive
        assert!((state.amplitudes[0].real - 1.0 / f64::sqrt(2.0)).abs() < 1e-10);
        
        // Second amplitude should become negative
        assert!((state.amplitudes[1].real + 1.0 / f64::sqrt(2.0)).abs() < 1e-10);
        
        // Both imaginary parts should be ~0
        assert!(state.amplitudes[0].imag.abs() < 1e-10);
        assert!(state.amplitudes[1].imag.abs() < 1e-10);
    }

    #[test]
    fn test_x_gate() {
        use crate::state::QuantumState;
        use crate::complex::Complex;
        
        // Test X flips |0⟩ to |1⟩
        let mut state = QuantumState::new(1);
        apply_x(&mut state, 0);
        
        // Should now be |1⟩
        assert!(state.amplitudes[0].magnitude_squared() < 1e-10); // |0⟩ amplitude is 0
        assert!((state.amplitudes[1].real - 1.0).abs() < 1e-10);  // |1⟩ amplitude is 1
        
        // Test X flips |1⟩ back to |0⟩
        apply_x(&mut state, 0);
        assert!((state.amplitudes[0].real - 1.0).abs() < 1e-10);  // Back to |0⟩
        assert!(state.amplitudes[1].magnitude_squared() < 1e-10);
    }

    #[test]
    fn test_y_gate() {
        use crate::state::QuantumState;
        
        // Test Y on |0⟩ gives i|1⟩
        let mut state = QuantumState::new(1);
        apply_y(&mut state, 0);
        
        // Should be i|1⟩, so amplitude[1] = 0 + 1i
        assert!(state.amplitudes[0].magnitude_squared() < 1e-10);
        assert!(state.amplitudes[1].real.abs() < 1e-10);
        assert!((state.amplitudes[1].imag - 1.0).abs() < 1e-10);
    }
}

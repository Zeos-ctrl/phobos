use crate::complex::Complex;

use rand::Rng;

/// The Quantum State Vector that holds the state of the system for 2^n qubits
#[derive(Clone)]
pub struct QuantumState {
    pub amplitudes: Vec<Complex>,
    pub num_qubits: usize,
}

impl QuantumState {
    /// Create state with n qubits, all in |0⟩
    ///
    /// # Arguments
    /// * `num_qubits` - Number of qubits to initalise the state with
    ///
    /// # Examples
    /// ```
    /// use phobos::QuantumState;
    ///
    /// let new_state = QuantumState::new(2);
    /// ```
    pub fn new(num_qubits: usize) -> QuantumState {
        let size = 2_usize.pow(num_qubits as u32);    
        let mut amplitudes = Vec::new();

        for i in 0..size {
            if i == 0 {
                amplitudes.push(Complex::new(1.0, 0.0));
            } else {
                amplitudes.push(Complex::new(0.0, 0.0));
            }
        }

        QuantumState { amplitudes, num_qubits }
    }

    /// Ensure sum of all amplitudes in the Quantum State equal 1, |amplitude|² = 1
    ///
    /// # Examples
    /// ```
    /// use phobos::QuantumState;
    ///
    /// let state = QuantumState::new(2);
    /// state.normalize();
    /// ```
    pub fn normalize(&mut self) {
        let mut sum = 0.0;
        
        // Calculate sum of |amplitude|²
        for amplitude in &self.amplitudes {
            sum += amplitude.magnitude_squared()
        }
        
        // Handle edge case
        if sum == 0.0 {
            return; // Can't normalize zero vector
        }
        
        let norm = f64::sqrt(sum);
        
        // Scale each amplitude
        for amplitude in &mut self.amplitudes {
            *amplitude = amplitude.scale(1.0 / norm);
        }
    }

    /// Measure all qubits, return bitstring (e.g., "01" for 2 qubits)
    ///
    /// # Examples
    /// ```
    /// use phobos::QuantumState;
    ///
    /// let state = QuantumState::new(2);
    /// let measurement = state.measure();
    /// ```
    pub fn measure(&self) -> String {
        let mut rng = rand::rng();
        let r: f64 = rng.random();
        
        let mut cumulative = 0.0;
        let mut chosen_index = 0;
        
        // Loop through amplitudes and find which index to measure
        for (index, amplitude) in self.amplitudes.iter().enumerate() {
            cumulative += amplitude.magnitude_squared();
            
            if r < cumulative {
                chosen_index = index;
                break;
            }
        }

        format!("{:0width$b}", chosen_index, width = self.num_qubits)
    }

}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_new_single_qubit() {
        let state = QuantumState::new(1);
        assert_eq!(state.num_qubits, 1);
        assert_eq!(state.amplitudes.len(), 2);
        // State should be |0⟩
        assert_eq!(state.amplitudes[0].real, 1.0);
        assert_eq!(state.amplitudes[0].imag, 0.0);
        assert_eq!(state.amplitudes[1].real, 0.0);
    }

    #[test]
    fn test_normalize() {
        let mut state = QuantumState::new(1);
        // Manually set unnormalized amplitudes
        state.amplitudes[0] = Complex::new(3.0, 0.0);
        state.amplitudes[1] = Complex::new(4.0, 0.0);
        
        state.normalize();
        
        // After normalization, |3|² + |4|² should become |0.6|² + |0.8|² = 1
        assert!((state.amplitudes[0].real - 0.6).abs() < 1e-10);
        assert!((state.amplitudes[1].real - 0.8).abs() < 1e-10);
    }

    #[test]
    fn test_measure_deterministic() {
        let state = QuantumState::new(1);
        // State is |0⟩, so measurement should always give "0"
        let result = state.measure();
        assert_eq!(result, "0");
        
        // Measure multiple times to be sure
        for _ in 0..10 {
            assert_eq!(state.measure(), "0");
        }
    }
}

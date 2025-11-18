use std::fmt;

/// Quantum gates supported by the simulator.
///
/// All gates are represented with their target qubit index (0-indexed).
/// Multi-qubit gates like CNOT also specify control qubits.
#[derive(Debug, Clone, Copy)]
pub enum Gate {
    /// Hadamard gate - creates superposition
    /// Transforms |0⟩ → (1/√2)(|0⟩ + |1⟩) and |1⟩ → (1/√2)(|0⟩ - |1⟩)
    Hadamard { target: usize },
    
    /// Controlled-NOT gate - flips target if control is |1⟩
    CNOT { control: usize, target: usize },
    
    /// Pauli-X gate - bit flip (quantum NOT)
    /// Transforms |0⟩ → |1⟩ and |1⟩ → |0⟩
    X { target: usize },
    
    /// Pauli-Y gate - bit flip with phase
    /// Transforms |0⟩ → i|1⟩ and |1⟩ → -i|0⟩
    Y { target: usize },
    
    /// Pauli-Z gate - phase flip
    /// Transforms |0⟩ → |0⟩ and |1⟩ → -|1⟩
    Z { target: usize },
    
    /// Identity gate - no operation (placeholder)
    I { target: usize },

    /// Controlled-Z gate - applies phase flip to target if control is |1⟩
    CZ { control: usize, target: usize },

    /// Measurement - collapses the target qubit to |0⟩ or |1⟩
    Measure { target: usize },

    /// Swap - Swaps the target qubits amplitudes
    Swap { qubit_a: usize, qubit_b: usize },

    /// Controlled phase rotation - applies a phase rotation to the target qubit only when
    /// both the control and target qubits are in the |1⟩ state.
    CPhase { control: usize, target: usize, angle: f64 }
}

/// A quantum circuit consisting of a sequence of gates applied to qubits.
///
/// Circuits are built by adding gates in sequence. When executed by a
/// [`crate::Simulator`], gates are applied in the order they were added.
pub struct Circuit {
    num_qubits: usize,
    gates: Vec<Gate>,
}

impl Circuit {
    /// Creates a new quantum circuit with the specified number of qubits.
    ///
    /// All qubits are initialized to |0⟩ when the circuit is executed.
    ///
    /// # Arguments
    /// * `num_qubits` - Number of qubits in the circuit
    ///
    /// # Examples
    /// ```
    /// use phobos::Circuit;
    /// let circuit = Circuit::new(2);
    /// ```
    pub fn new(num_qubits: usize) -> Self {
        Circuit { num_qubits, gates: Vec::new() }
    }
    
    /// Adds a gate to the circuit.
    ///
    /// Gates are executed in the order they are added.
    ///
    /// # Arguments
    /// * `gate` - The gate to add to the circuit
    ///
    /// # Examples
    /// ```
    /// use phobos::{Circuit, Gate};
    /// let mut circuit = Circuit::new(2);
    /// circuit.add_gate(Gate::Hadamard { target: 0 });
    /// circuit.add_gate(Gate::CNOT { control: 0, target: 1 });
    /// ```
    pub fn add_gate(&mut self, gate: Gate) {
        self.gates.push(gate)
    }

    /// Applies the Quantum Fourier Transform to a range of qubits.
    ///
    /// The QFT decomposes quantum states based on their phase rotation patterns,
    /// analogous to how a classical FFT decomposes signals into frequency components.
    ///
    /// # Arguments
    /// * `qubit_range` - The range of qubits to apply QFT to (e.g., 0..3 for qubits 0,1,2)
    ///
    /// # Examples
    /// ```
    /// use phobos::Circuit;
    ///
    /// let mut circuit = Circuit::new(4);
    /// circuit.apply_qft(0..4);  // Apply QFT to all 4 qubits
    /// ```
    pub fn apply_qft(&mut self, qubit_range: std::ops::Range<usize>) {
        let qubits: Vec<usize> = qubit_range.collect();
        let n = qubits.len();
        
        // Process each qubit
        for i in 0..n {
            let current_qubit = qubits[i];
            
            self.add_gate(Gate::Hadamard { target: current_qubit });
            
            // Apply controlled rotations from all qubits after current_qubit
            for j in i + 1..n {
                let control_qubit = qubits[j];
                let target_qubit = qubits[i];
                
                // compute angle and add CPhase gate
                let angle = 2.0 * std::f64::consts::PI / (1_u32 << (j - i + 1)) as f64;
                self.add_gate(Gate::CPhase{ control: control_qubit, target: target_qubit, angle: angle})
            }
        }
        
        // Swap qubits to reverse order
        for i in 0..n/2 {
            self.add_gate(Gate::Swap { qubit_a: qubits[i], qubit_b: qubits[n - 1 -i] });
        }
    }
    
    /// Returns the number of qubits in the circuit.
    pub fn num_qubits(&self) -> usize {
        self.num_qubits
    }
    
    /// Returns a slice of all gates in the circuit.
    pub fn gates(&self) -> &[Gate] {
        &self.gates
    }
}

impl fmt::Display for Circuit {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        // Create a vector of strings, one per qubit
        let mut lines: Vec<String> = Vec::new();
        
        // Initialize each line with the qubit label
        for i in 0..self.num_qubits {
            lines.push(format!("q{}: ─", i));
        }
        
        // Process each gate
        for gate in self.gates() {
            match gate {
                Gate::Hadamard { target } => {
                    for i in 0..self.num_qubits {
                        match i {
                            i if i == *target => {
                                lines[i].push('H');
                                lines[i].push('─');
                            },
                            _ => {
                                lines[i].push('─');
                                lines[i].push('─');
                            }
                        }
                    }
                },
                Gate::CNOT { control, target } => {
                    for i in 0..self.num_qubits {
                        match i {
                            i if i == *target => {
                                lines[i].push('C');
                                lines[i].push('─')
                            },
                            i if i == *control => {
                                lines[i].push('●');
                                lines[i].push('─')
                            },
                            _ => {
                                lines[i].push('─');
                                lines[i].push('─')
                            }
                        }
                    }
                },
                Gate::X { target } => {
                    for i in 0..self.num_qubits {
                        match i {
                            i if i == *target => {
                                lines[i].push('X');
                                lines[i].push('─');
                            },
                            _ => {
                                lines[i].push('─');
                                lines[i].push('─');
                            }
                        }
                    }
                },
                Gate::Y { target } => {
                    for i in 0..self.num_qubits {
                        match i {
                            i if i == *target => {
                                lines[i].push('Y');
                                lines[i].push('─');
                            },
                            _ => {
                                lines[i].push('─');
                                lines[i].push('─');
                            }
                        }
                    }
                },
                Gate::Z { target } => {
                    for i in 0..self.num_qubits {
                        match i {
                            i if i == *target => {
                                lines[i].push('Z');
                                lines[i].push('─');
                            },
                            _ => {
                                lines[i].push('─');
                                lines[i].push('─');
                            }
                        }
                    }
                },
                Gate::I { target } => {
                    for i in 0..self.num_qubits {
                        match i {
                            i if i == *target => {
                                lines[i].push('I');
                                lines[i].push('─');
                            },
                            _ => {
                                lines[i].push('─');
                                lines[i].push('─');
                            }
                        }
                    }
                },
                Gate::CZ { control, target } => {
                    for i in 0..self.num_qubits {
                        match i {
                            i if i == *target => {
                                lines[i].push('Z');
                                lines[i].push('─');
                            },
                            i if i == *control => {
                                lines[i].push('●');
                                lines[i].push('─')
                            },
                            _ => {
                                lines[i].push('─');
                                lines[i].push('─');
                            }
                        }
                    }
                },
                Gate::Swap { qubit_a, qubit_b } => {
                    for i in 0..self.num_qubits {
                        match i {
                            i if i == *qubit_a => {
                                lines[i].push('S');
                                lines[i].push('─');
                            },
                            i if i == *qubit_b => {
                                lines[i].push('S');
                                lines[i].push('─')
                            },
                            _ => {
                                lines[i].push('─');
                                lines[i].push('─');
                            }
                        }
                    }
                },
                Gate::CPhase { control, target, .. } => {
                    for i in 0..self.num_qubits {
                        match i {
                            i if i == *target => {
                                lines[i].push('P');
                                lines[i].push('─');
                            },
                            i if i == *control => {
                                lines[i].push('●');
                                lines[i].push('─')
                            },
                            _ => {
                                lines[i].push('─');
                                lines[i].push('─');
                            }
                        }
                    }
                },
                Gate::Measure { target } => {
                    for i in 0..self.num_qubits {
                        match i {
                            i if i == *target => {
                                lines[i].push('M');
                                lines[i].push('─');
                            },
                            _ => {
                                lines[i].push('─');
                                lines[i].push('─');
                            }
                        }
                    }
                },
            }
        }
        
        for line in lines {
            writeln!(f, "{}", line)?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_new_circuit() {
        let circuit = Circuit::new(2);
        assert_eq!(circuit.num_qubits(), 2);
        assert_eq!(circuit.gates().len(), 0);
    }

    #[test]
    fn test_adding_gates() {
        let mut circuit = Circuit::new(5);
        circuit.add_gate(Gate::Hadamard { target: 0 });
        circuit.add_gate(Gate::CNOT { control: 0, target: 1 });
        circuit.add_gate(Gate::X { target: 0 });
        circuit.add_gate(Gate::Y { target: 0 });
        circuit.add_gate(Gate::Z { target: 0 });
        
        assert_eq!(circuit.gates().len(), 5);
        
        // Check first gate
        match circuit.gates()[0] {
            Gate::Hadamard { target } => assert_eq!(target, 0),
            _ => panic!("Expected Hadamard gate"),
        }
        
        // Check second gate
        match circuit.gates()[1] {
            Gate::CNOT { control, target } => {
                assert_eq!(control, 0);
                assert_eq!(target, 1);
            }
            _ => panic!("Expected CNOT gate"),
        }

        // Check third gate
        match circuit.gates()[2] {
            Gate::X { target } => assert_eq!(target, 0),
            _ => panic!("Expected X gate"),
        }

        // Check fourth gate
        match circuit.gates()[3] {
            Gate::Y { target } => assert_eq!(target, 0),
            _ => panic!("Expected Y gate"),
        }

        // Check fifth gate
        match circuit.gates()[4] {
            Gate::Z { target } => assert_eq!(target, 0),
            _ => panic!("Expected Z gate"),
        }
    }

   #[test]
    fn test_qft_circuit() {
        let mut circuit = Circuit::new(4);
        circuit.apply_qft(0..4);
        
        println!("QFT Circuit:");
        println!("{}", circuit);
        
        assert_eq!(circuit.gates().len(), 12);
    } 
}

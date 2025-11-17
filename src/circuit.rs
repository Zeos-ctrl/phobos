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
}

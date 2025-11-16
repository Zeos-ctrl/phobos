use std::fmt;

#[derive(Debug, Clone, Copy)]
pub enum Gate {
    Hadamard { target: usize },
    CNOT { control: usize, target: usize },
}

pub struct Circuit {
    num_qubits: usize,
    gates: Vec<Gate>,
}

impl Circuit {
    pub fn new(num_qubits: usize) -> Self {
        Circuit { num_qubits, gates: Vec::new() }
    }

    pub fn add_gate(&mut self, gate: Gate) {
        self.gates.push(gate)
    }

    pub fn num_qubits(&self) -> usize {
        self.num_qubits
    }

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
                }
                Gate::CNOT { control, target } => {
                    for i in 0..self.num_qubits {
                        match i {
                            i if i == *target => {
                                lines[i].push('X');
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
                }
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
        let mut circuit = Circuit::new(2);
        circuit.add_gate(Gate::Hadamard { target: 0 });
        circuit.add_gate(Gate::CNOT { control: 0, target: 1 });
        
        assert_eq!(circuit.gates().len(), 2);
        
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
    }
}

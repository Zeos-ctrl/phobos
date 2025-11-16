use crate::circuit::{Circuit, Gate};
use crate::state::QuantumState;
use crate::Complex;
use crate::gates;

use std::fmt;

pub struct Simulator;

pub struct TraceStep {
    pub description: String,
    pub state: QuantumState,
}

pub struct ExecutionTrace {
    pub steps: Vec<TraceStep>,
}

impl fmt::Display for ExecutionTrace {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        for (step_num, step) in self.steps.iter().enumerate() {
            writeln!(f, "Step {}: {}", step_num, step.description)?;
            
            // Convert state to ket notation
            let ket = state_to_ket(&step.state);
            writeln!(f, "  State: {}", ket)?;
            writeln!(f)?; // blank line between steps
        }
        Ok(())
    }
}

fn state_to_ket(state: &QuantumState) -> String {
    let mut terms = Vec::new();
    
    for (index, amplitude) in state.amplitudes.iter().enumerate() {
        // Skip near-zero amplitudes
        if amplitude.magnitude_squared() < 1e-10 {
            continue;
        }
        
        // Convert index to binary basis state
        let basis = format!("{:0width$b}", index, width = state.num_qubits);
        
        // Format the coefficient (just the real part for now, assuming imaginary ≈ 0)
        let coeff = amplitude.real;
        
        terms.push(format!("{:.4}|{}⟩", coeff, basis));
    }
    
    terms.join(" + ")
}

impl Simulator {
    pub fn new() -> Self {
        Simulator
    }

    // Run circuit n_shots times, return measurement results
    pub fn run(&self, circuit: &Circuit, n_shots: usize) -> Vec<String> {
        let mut results = Vec::new();

        for _ in 0..n_shots {
            results.push(self.run_once(&circuit));
        }

        results
    }

    // Helper: Execute circuit once and return final measurement
    fn run_once(&self, circuit: &Circuit) -> String {
        let mut state = QuantumState::new(circuit.num_qubits());
        for gate in circuit.gates() {
            match gate {
                Gate::Hadamard { target } => {
                    gates::apply_hadamard(&mut state, *target);
                },
                Gate::CNOT { control, target } => {
                    gates::apply_cnot(&mut state, *control, *target);
                },
                Gate::X { target } => {
                    gates::apply_x(&mut state, *target);
                },
                Gate::Y { target } => {
                    gates::apply_y(&mut state, *target);
                },
                Gate::Z { target } => {
                    gates::apply_z(&mut state, *target);
                },
                Gate::I { target } => {
                    gates::apply_identity(&mut state, *target);
                }
            }
        }
        state.measure()
    }

    pub fn run_with_trace(&self, circuit: &Circuit) -> ExecutionTrace {
        let mut trace = ExecutionTrace { steps: Vec::new() };
        let mut state = QuantumState::new(circuit.num_qubits());
        
        // Capture initial state
        trace.steps.push(TraceStep {
            description: "Initial state".to_string(),
            state: state.clone(), // We need to clone the state
        });
        
        // Apply each gate and capture state after each
        for gate in circuit.gates() {
            match gate {
                Gate::Hadamard { target } => {
                    gates::apply_hadamard(&mut state, *target);
                    trace.steps.push(TraceStep {
                        description: format!("Applied Hadamard to qubit {}", *target),
                        state: state.clone()
                    });
                },
                Gate::CNOT { control, target } => {
                    gates::apply_cnot(&mut state, *control, *target);
                    trace.steps.push(TraceStep {
                        description: format!("Applied CNOT (control: {}, target: {})", *control, *target),
                        state: state.clone()
                    });
                },
                Gate::X { target } => {
                    gates::apply_x(&mut state, *target);
                    trace.steps.push(TraceStep {
                        description: format!("Applied X gate to qubit {}", *target),
                        state: state.clone()
                    });
                },
                Gate::Y { target } => {
                    gates::apply_y(&mut state, *target);
                    trace.steps.push(TraceStep {
                        description: format!("Applied Y gate to qubit {}", *target),
                        state: state.clone()
                    });
                },
                Gate::Z { target } => {
                    gates::apply_z(&mut state, *target);
                    trace.steps.push(TraceStep {
                        description: format!("Applied Z gate to qubit {}", *target),
                        state: state.clone()
                    });
                },
                Gate::I { target } => {
                    gates::apply_identity(&mut state, *target);
                    trace.steps.push(TraceStep {
                        description: format!("Applied I gate to qubit {}", *target),
                        state: state.clone()
                    });
                },
            }
        }
        let result = state.measure();

        let measured_index = usize::from_str_radix(&result, 2).unwrap();

        // Set all amplitudes to 0, then set measured state to 1
        for i in 0..state.amplitudes.len() {
            if i == measured_index {
                state.amplitudes[i] = Complex::new(1.0, 0.0);
            } else {
                state.amplitudes[i] = Complex::new(0.0, 0.0);
            }
        }

        trace.steps.push(TraceStep {
            description: format!("Measurement result: {}", result),
            state: state.clone(),
        });

        trace
    }
}

use crate::circuit::{Circuit, Gate};
use crate::state::QuantumState;
use crate::Complex;
use crate::gates;

use std::fmt;

/// Simulator for the quantum circuit
pub struct Simulator;

impl Default for Simulator {
    fn default() -> Self {
        Simulator::new()
    }
}

/// Struct to hold information of the current step for tracing.
/// Holds a description of the step as well as the Quantum State.
pub struct TraceStep {
    pub description: String,
    pub state: QuantumState,
}

/// Holds multiple TraceSteps to be iterated over when read.
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
    /// Returns a new Simulator object to run a circuit.
    ///
    /// # Examples
    /// ```
    /// use phobos::Simulator;
    ///
    /// let sim = Simulator::new();
    /// ```
    pub fn new() -> Self {
        Simulator
    }

    /// Run circuit n_shots times, return measurement results.
    ///
    /// # Arguments
    /// * `&Circuit` - Reference to a circuit to run
    /// * `n_shots` - Number of times to run the circuit
    ///
    /// # Examples
    /// ```
    /// use phobos::Simulator;
    /// use phobos::Circuit;
    ///
    /// let circuit = Circuit::new(1);
    /// let sim = Simulator::new();
    ///
    /// let results = sim.run(&circuit, 100);
    /// ```
    pub fn run(&self, circuit: &Circuit, n_shots: usize) -> Vec<String> {
        let mut results = Vec::new();

        for _ in 0..n_shots {
            results.push(self.run_once(circuit));
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
                Gate::CZ { control, target } => {
                    gates::apply_cz(&mut state, *control, *target);
                },
                Gate::Swap { qubit_a, qubit_b } => {
                    gates::apply_swap(&mut state, *qubit_a, *qubit_b);
                },
                Gate::CPhase { control, target, angle } => {
                    gates::apply_cphase(&mut state, *control, *target, *angle);
                },
                Gate::CSwap { control, qubit_a, qubit_b } => {
                    gates::apply_cswap(&mut state, *control, *qubit_a, *qubit_b);
                },
                Gate::I { target } => {
                    gates::apply_identity(&mut state, *target);
                },
                Gate::U { target, matrix } => {
                    gates::apply_u(&mut state, *target, matrix);
                },
                Gate::Measure { target } => {
                    gates::measure_qubit(&mut state, *target);
                },
            }
        }
        state.measure()
    }

    /// Run circuit with the trace active to see the journey of the quantum state through
    /// the circuit.
    ///
    /// # Arguments
    /// * `&Circuit` - Reference to a circuit to run
    ///
    /// # Examples
    /// ```
    /// use phobos::Simulator;
    /// use phobos::Circuit;
    ///
    /// let circuit = Circuit::new(1);
    /// let sim = Simulator::new();
    ///
    /// let trace = sim.run_with_trace(&circuit);
    /// println!("Circuit Trace:");
    /// println!("{}", trace);
    /// ```
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
                Gate::CZ { control, target } => {
                    gates::apply_cz(&mut state, *control, *target);
                    trace.steps.push(TraceStep {
                        description: format!("Applied CZ gate (control: {}, target: {})", *control, *target),
                        state: state.clone()
                    });
                },
                Gate::Swap { qubit_a, qubit_b } => {
                    gates::apply_swap(&mut state, *qubit_a, *qubit_b);
                    trace.steps.push(TraceStep {
                        description: format!("Applied Swap gate (A: {}, B: {})", *qubit_a, *qubit_b),
                        state: state.clone()
                    });
                },
                Gate::CPhase { control, target, angle } => {
                    gates::apply_cphase(&mut state, *control, *target, *angle);
                    trace.steps.push(TraceStep {
                        description: format!("Applied Controlled Phase Rotation (control: {}, target: {}, angle: {})", *control, *target, *angle),
                        state: state.clone()
                    });
                },
                Gate::CSwap { control, qubit_a, qubit_b } => {
                    gates::apply_cswap(&mut state, *control, *qubit_a, *qubit_b);
                    trace.steps.push(TraceStep {
                        description: format!("Applied Controlled Swap (control: {}, A: {}, B: {})", *control, *qubit_a, *qubit_b),
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
                Gate::U { target, matrix } => {
                    gates::apply_u(&mut state, *target, matrix);
                    trace.steps.push(TraceStep {
                        description: format!("Applied U gate to qubit {}, with matrix {:?}", *target, matrix),
                        state: state.clone()
                    });
                },
                Gate::Measure { target } => {
                    let result = gates::measure_qubit(&mut state, *target);
                    trace.steps.push(TraceStep {
                        description: format!("Measured qubit {}: result = {}", target, result),
                        state: state.clone(),
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

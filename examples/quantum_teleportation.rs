use phobos::{Circuit, Gate, Simulator, plot_histogram_qubit};

fn quantum_teleportation(qubit_to_send_op: &str, num_shots: usize) {
    // Create 3-qubit circuit
    // Q0: Alice's state to teleport
    // Q1: Alice's control qubit
    // Q2: Bob's control qubit
    let mut circuit = Circuit::new(3);
    
    // Prepare Q0 with the state to teleport
    match qubit_to_send_op {
        "H" => circuit.add_gate(Gate::Hadamard { target: 0 }),
        "X" => circuit.add_gate(Gate::X { target: 0 }),
        "Y" => circuit.add_gate(Gate::Y { target: 0 }),
        "Z" => circuit.add_gate(Gate::Z { target: 0 }),
        "I" => circuit.add_gate(Gate::I { target: 0 }),
        _ => panic!("Unsupported gate: {}", qubit_to_send_op),
    }
    
    // Create Bell state between Q1 and Q2
    circuit.add_gate(Gate::Hadamard { target: 1 });
    circuit.add_gate(Gate::CNOT{ control: 1, target: 2 });
    
    // CNOT with Q0 as control, Q1 as target
    circuit.add_gate(Gate::CNOT{ control: 0, target: 1 });
    
    // Hadamard on Q0, then measure Q0 and Q1
    circuit.add_gate(Gate::Hadamard { target: 0 });
    circuit.add_gate(Gate::Measure { target: 0 });
    circuit.add_gate(Gate::Measure { target: 1 });
    
    // Apply corrections to Q2 based on measurements
    circuit.add_gate(Gate::CNOT { control: 1, target: 2 });
    circuit.add_gate(Gate::CZ {control: 0, target: 2});
    
    // Print circuit
    println!("Circuit:");
    println!("{}", circuit);
    
    // Run simulation
    let sim = Simulator::new();
    let results = sim.run(&circuit, num_shots);
    
    // Show histogram
    println!("\nMeasurement results ({} shots):", num_shots);
    plot_histogram_qubit(&results, circuit.num_qubits(), Some(2));
}

fn main() {
    quantum_teleportation("H", 1000);
}

use phobos::{Circuit, Gate, Simulator, plot_histogram};

fn main() {
    // Create circuit with 3 qubits
    let mut circuit = Circuit::new(3);
    
    // Prepare |101⟩ state
    println!("Prepareing |101⟩ state");
    circuit.add_gate(Gate::X { target: 0 });
    circuit.add_gate(Gate::X { target: 2 });
    
    println!("=== Initial State Preparation ===");
    println!("{}", circuit);
    
    // Apply QFT
    circuit.apply_qft(0..3);
    println!("\n=== After QFT ===");
    println!("{}", circuit);
    
    // Apply inverse QFT
    circuit.apply_inverse_qft(0..3);
    println!("\n=== After Inverse QFT ===");
    println!("{}", circuit);
    
    // Run simulation
    let simulator = Simulator::new();
    let results = simulator.run(&circuit, 1000);
    
    println!("\n=== Measurement Results (1000 shots) ===");
    // Count outcomes
    let mut count_101 = 0;
    let mut count_0 = 0;
    
    for result in &results {
        match result.as_str() {
            "101" => count_101 += 1,
            _ => { count_0 += 1}
        }
    }

    println!("\nResults from 1000 measurements:");
    println!("101: {} ({:.1}%)", count_101, count_101 as f64 / 10.0);
    println!("Others: {} ({:.1}%)", count_0, count_0 as f64 / 10.0);

    println!("\nProbability Distribution:");
    plot_histogram(&results, circuit.num_qubits());
}

use phobos::{Circuit, Gate, Simulator, QPEUnitary, plot_histogram};

fn main() {
    let mut circuit = Circuit::new(3);  // 2 readout + 1 target
    
    // Prepare eigenstate |1⟩ in target qubit
    circuit.add_gate(Gate::X { target: 2 });
    
    println!("=== QPE: Estimating phase of Z gate acting on |1⟩ ===\n");
    
    // Apply QPE
    circuit.apply_qpe(0..2, 2, QPEUnitary::Z);
    
    println!("Circuit:");
    println!("{}", circuit);
    
    // Run simulation
    let simulator = Simulator::new();
    let results = simulator.run(&circuit, 1000);
    
    println!("\n=== Measurement Results (1000 shots) ===");

    // Count outcomes
    let mut count_101 = 0;
    let mut count_other = 0;

    for result in &results {
        match result.as_str() {
            "101" => count_101 += 1,
            _ => count_other += 1,
        }
    }

    println!("\nResults from 1000 measurements:");
    println!("101 (readout=10, target=1): {} ({:.1}%)", count_101, count_101 as f64 / 10.0);
    println!("Other: {} ({:.1}%)", count_other, count_other as f64 / 10.0);

    println!("\nProbability Distribution:");
    plot_histogram(&results, circuit.num_qubits());

    println!("\nPhase Estimation Result:");
    println!("Readout qubits measured: |10⟩ = 0.10 in binary = 0.5 in decimal");
}

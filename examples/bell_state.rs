use phobos::{Circuit, Gate, Simulator, plot_histogram};

fn main() {
    // Create a Bell state circuit
    let mut circuit = Circuit::new(2);
    circuit.add_gate(Gate::Hadamard { target: 0 });
    circuit.add_gate(Gate::CNOT { control: 0, target: 1 });
    
    println!("Circuit:");
    println!("{}", circuit);
    
    // Run the circuit
    let sim = Simulator::new();
    let results = sim.run(&circuit, 1000);
    
    // Count outcomes
    let mut count_00 = 0;
    let mut count_11 = 0;
    
    for result in &results {
        match result.as_str() {
            "00" => count_00 += 1,
            "11" => count_11 += 1,
            _ => {}
        }
    }

    println!("\nResults from 1000 measurements:");
    println!("00: {} ({:.1}%)", count_00, count_00 as f64 / 10.0);
    println!("11: {} ({:.1}%)", count_11, count_11 as f64 / 10.0);

    println!("\nProbability Distribution:");
    plot_histogram(&results, circuit.num_qubits());
}

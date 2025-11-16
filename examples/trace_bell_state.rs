use phobos::{Circuit, Gate, Simulator};

fn main() {
    let mut circuit = Circuit::new(2);
    circuit.add_gate(Gate::Hadamard { target: 0 });
    circuit.add_gate(Gate::CNOT { control: 0, target: 1 });
    
    let sim = Simulator::new();
    let trace = sim.run_with_trace(&circuit);

    println!("Circuit:");
    println!("{}", circuit);
    
    println!("Circuit Trace:");
    println!("{}", trace);
}

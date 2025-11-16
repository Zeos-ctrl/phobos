use phobos::{Circuit, Gate};

fn main() {
    // Create a simple Bell state circuit
    let mut circuit = Circuit::new(2);
    circuit.add_gate(Gate::Hadamard { target: 0 });
    circuit.add_gate(Gate::CNOT { control: 0, target: 1 });
    
    println!("{}", circuit);
}

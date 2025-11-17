# Phobos - Quantum Circuit Simulator
A from-scratch quantum circuit simulator built in Rust for learning and experimentation with quantum computing.

## Installation

Add this to your Cargo.toml:
``` toml
[dependencies]
phobos = "0.1.0"
```

Or clone and build from source:

``` bash
git clone https://github.com/yourusername/phobos.git
cd phobos
cargo build --release
```

## Quick Start

Here's a simple example creating a Bell state (maximally entangled pair):

``` rust
use phobos::{Circuit, Gate, Simulator};

fn main() {
    // Create a 2-qubit circuit
    let mut circuit = Circuit::new(2);
    
    // Add gates
    circuit.add_gate(Gate::Hadamard { target: 0 });
    circuit.add_gate(Gate::CNOT { control: 0, target: 1 });
    
    // Run the circuit 1000 times
    let sim = Simulator::new();
    let results = sim.run(&circuit, 1000);
    
    // Results will be roughly 50% "00" and 50% "11"
    println!("Measurement results: {:?}", results);
}
```

## Visualising Circuits

Phobos includes some visualising functions. We can print the circuit diagram with:

``` rust
println!("{}", circuit);
```

Output:
``` bash
Circuit:
q0: ─H─●─
q1: ───C─
```

We can see there is a Hadamard gate on qubit 0, and then we use it as a target for a CNOT gate on qubit 1.

With phobos you can trace the state of the state vector through the circuit:

``` rust
use phobos::{Circuit, Gate, Simulator};

fn main() {
    let mut circuit = Circuit::new(2);
    circuit.add_gate(Gate::Hadamard { target: 0 });
    circuit.add_gate(Gate::CNOT { control: 0, target: 1 });
    
    let sim = Simulator::new();
    let trace = sim.run_with_trace(&circuit);

    println!("Circuit Trace:");
    println!("{}", trace);
}
```

Output:
``` bash
Circuit Trace:
Step 0: Initial state
  State: 1.0000|00⟩

Step 1: Applied Hadamard to qubit 0
  State: 0.7071|00⟩ + 0.7071|01⟩

Step 2: Applied CNOT (control: 0, target: 1)
  State: 0.7071|00⟩ + 0.7071|11⟩

Step 3: Measurement result: 00
  State: 1.0000|00⟩
```

We can also visualise the histogram of the output of multiple measurements:

``` rust
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
```

Output:
```bash
Results from 1000 measurements:
00: 499 (49.9%)
11: 501 (50.1%)

Probability Distribution:
00: ❚❚❚❚❚❚❚❚❚❚❚❚❚❚❚❚❚❚❚❚❚❚❚❚❚❚❚❚❚❚❚❚❚❚❚❚❚❚❚❚❚❚❚❚❚❚❚❚❚
01:
10:
11: ❚❚❚❚❚❚❚❚❚❚❚❚❚❚❚❚❚❚❚❚❚❚❚❚❚❚❚❚❚❚❚❚❚❚❚❚❚❚❚❚❚❚❚❚❚❚❚❚❚❚
```

## Examples

You can find example circuits in `Examples/Circuit_name.rs`.

use phobos::{Circuit, Gate, Simulator};
use std::collections::HashMap;

fn gcd(mut a: i32, mut b: i32) -> i32 {
    while b != 0 {
        let temp = b;
        b = a % b;
        a = temp;
    }
    a
}

fn periodic_oracle(circuit: &mut Circuit, 
                               output_qubits: &[usize], 
                               input_qubits: &[usize]) {
    // For each output qubit, apply the oracle U^(2^k) times
    for (k, &output_qubit) in output_qubits.iter().enumerate() {
        let power = 2_u32.pow((output_qubits.len() - 1 - k) as u32);
        
        // Apply the oracle 'power' times
        for _ in 0..power {
            circuit.add_gate(Gate::CSwap { control: output_qubit, qubit_a: input_qubits[2], qubit_b: input_qubits[3] });
            circuit.add_gate(Gate::CSwap { control: output_qubit, qubit_a: input_qubits[1], qubit_b: input_qubits[2] });
            circuit.add_gate(Gate::CSwap { control: output_qubit, qubit_a: input_qubits[0], qubit_b: input_qubits[1] });

            for j in 0..4 {
                circuit.add_gate(Gate::CNOT { control: output_qubit, target: input_qubits[j]})
            }
        }
    }
}

fn find_period_from_phases(phase_counts: &HashMap<String, usize>, num_output_qubits: usize) -> Option<u32> {
    let mut best_phase = 0.0;
    let mut best_s = 0;
    let mut best_r = 1;
    
    println!("\nAnalyzing all measurements:");
    println!("State    | Phase   | Rational");
    println!("---------|---------|----------");
    
    // Loop through ALL measurements
    for (bits, _count) in phase_counts.iter() {
        // Convert bits to phase
        let value = usize::from_str_radix(bits, 2).unwrap();
        let phase = value as f64 / (1 << num_output_qubits) as f64;
        
        // Find rational approximation s/r with denominator <= 15
        let denom_limit = 15;
        let mut this_best_r = 1;
        let mut this_best_s = 0;
        let mut this_best_error = f64::MAX;
        
        for r in 1..=denom_limit {
            let s = (phase * r as f64).round() as u32;
            let approx = s as f64 / r as f64;
            let error = (phase - approx).abs();
            
            if error < this_best_error {
                this_best_error = error;
                this_best_r = r;
                this_best_s = s;
            }
        }
        
        println!("{:8} | {:.4}  | {}/{}", bits, phase, this_best_s, this_best_r);
        
        // Keep track of maximum phase
        if phase > best_phase {
            best_phase = phase;
            best_s = this_best_s;
            best_r = this_best_r;
        }
    }
    
    println!("\nMaximum phase: {:.4} = {}/{}", best_phase, best_s, best_r);
    
    // Check if fraction is in form (r-1)/r
    if best_r - best_s == 1 {
        println!("Period found: r = {}", best_r);
        Some(best_r)
    } else {
        println!("Cannot determine period from this measurement (fraction not in (r-1)/r form)");
        None
    }
}

fn main() {
    println!("Shor's Period Finding: a=7, N=15");
    let a: i32 = 7;
    let n: i32 = 15;
    
    let num_output_qubits = 4;  // Precision bits
    let num_input_qubits = 4;   // For domain size 16
    let total_qubits = num_output_qubits + num_input_qubits;
    
    let mut circuit = Circuit::new(total_qubits);
    
    // Hadamards on output qubits (0-3)
    for i in 0..4 {
        circuit.add_gate( Gate::Hadamard { target: i });
    }
    
    // Prepare input qubits as |0001⟩ (qubit 7 = 1)
    circuit.add_gate( Gate::X { target: total_qubits - 1});
    
    // Apply periodic oracle (this is the new part!)
    periodic_oracle(&mut circuit, &[0, 1, 2, 3], &[4, 5, 6, 7]);
    
    // Inverse QFT on output qubits
    circuit.apply_inverse_qft(0..4);
    
    // Run simulation
    let simulator = Simulator::new();
    let results = simulator.run(&circuit, 1000);

    let mut phase_counts: std::collections::HashMap<String, usize> = std::collections::HashMap::new();

    for result in &results {
        // Take first 4 characters (output qubits)
        let output_bits = &result[..4];
        *phase_counts.entry(output_bits.to_string()).or_insert(0) += 1;
    }

    let mut candidate_r = find_period_from_phases(&phase_counts, num_output_qubits).unwrap();
    println!("Measured period: r = {}", candidate_r);

    while candidate_r % 2 == 0 {
        let half_r = candidate_r / 2;
        let a_to_half_r = a.pow(half_r);
        
        let factor1 = gcd(a_to_half_r - 1, n);
        let factor2 = gcd(a_to_half_r + 1, n);
        
        println!("\nFactorization:");
        println!("a^(r/2) = {}^{} = {}", a, half_r, a_to_half_r);
        println!("gcd({} - 1, {}) = gcd({}, {}) = {}", a_to_half_r, n, a_to_half_r - 1, n, factor1);
        println!("gcd({} + 1, {}) = gcd({}, {}) = {}", a_to_half_r, n, a_to_half_r + 1, n, factor2);
        

        if (factor1 * factor2 == n) && (1 < factor1 && factor1 < n) && (1 < factor2 && factor2 < n) {
            println!("factor 1 * factor 2 = {} ({} * {} = {})", n, factor1, factor2, n);
            break;
        } else {
            println!("Trying r = {} / 2 = {}", candidate_r, candidate_r / 2);
            candidate_r = candidate_r / 2;
        }
    }
}

use std::collections::HashMap;

/// Utility function to create a histogram of the results of a quantum circuit
/// simulation
///
/// # Arguments
/// * `results` - Result vector from the simulator
/// * `num_qubits` - Number of qubits in the system to get the depth of the histogram
///
/// # Examples
/// ```
/// use phobos::Circuit;
/// use phobos::Simulator;
/// use phobos::plot_histogram;
///
/// let circuit = Circuit::new(1);
/// let sim = Simulator::new();
/// let results = sim.run(&circuit, 1000);
///
/// println!("\nProbability Distribution:");
/// plot_histogram(&results, circuit.num_qubits()); 
/// ```
pub fn plot_histogram(results: &[String], num_qubits: usize) {
    let mut frequency = HashMap::new();

    // Get the state value as a key and iter over vector
    for state in results.iter() {
        *frequency.entry(state).or_insert(0) += 1;
    }

    let total = results.len() as f64;
    let num_states = 2_usize.pow(num_qubits as u32);

    for i in 0..num_states {
        // Convert i to binary string
        let binary = format!("{:0width$b}", i, width = num_qubits);
        // Look up count in frequency HashMap
        let count = frequency.get(&binary).unwrap_or(&0);
        // Calculate percentage
        let percentage = (*count as f64 / total) * 100.0;
        // Build bar string
        let bar = "❚".repeat(percentage.floor() as usize);
        println!("{}: {}", binary, bar);
    }
}

/// Plots a histogram for a specific qubit's measurement results.
///
/// Displays a horizontal bar chart showing the distribution of a single qubit's outcomes.
///
/// # Arguments
/// * `results` - Vector of measurement results (binary strings)
/// * `num_qubits` - Total number of qubits in the system
/// * `target_qubit` - Index of the qubit to display (0-indexed)
pub fn plot_histogram_qubit(results: &[String], num_qubits: usize, target_qubit: Option<usize>) {
    let mut frequency = HashMap::new();
    
    for state in results.iter() {
        let key = if let Some(qubit) = target_qubit {
            let bit_index = num_qubits - 1 - qubit;
            state.chars().nth(bit_index).unwrap().to_string()
        } else {
            state.clone()
        };
        *frequency.entry(key).or_insert(0) += 1;
    }
    
    let total = results.len() as f64;
    let (num_states, width) = if target_qubit.is_some() {
        (2, 1)
    } else {
        (2_usize.pow(num_qubits as u32), num_qubits)
    };
    
    for i in 0..num_states {
        let binary = format!("{:0width$b}", i, width = width);
        let count = frequency.get(&binary).unwrap_or(&0);
        let percentage = (*count as f64 / total) * 100.0;
        let bar = "❚".repeat(percentage.floor() as usize);
        println!("{}: {}", binary, bar);
    }
}

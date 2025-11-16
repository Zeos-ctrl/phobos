use std::collections::HashMap;

pub fn plot_histogram(results: &Vec<String>, num_qubits: usize) {
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

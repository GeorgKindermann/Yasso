mod yasso20;

use nalgebra::SVector;
use std::time::Instant;

fn main() {
    let raw_theta = vec![0.51,5.19,0.13,0.1,0.5,0.,1.,1.,0.99,0.,0.,0.,0.,0.,0.163,0.,-0.,0.,0.,0.,0.,0.158,-0.002,0.17,-0.005,0.067,-0.,-1.44,-2.0,-6.9,0.0042,0.0015,-2.55,1.24,0.25];
    let theta = yasso20::get_theta(raw_theta);

    let mut init = SVector::<f64, 5>::zeros();
    // Nutze from_row_slice für die Erstellung
    let infall = SVector::<f64, 5>::from_row_slice(&[0.5, 0.1, 0.1, 0.2, 0.0]);
    let avg_t = SVector::<f64, 12>::from_row_slice(&[-2., 2., 5., 10., 15., 20., 22., 21., 16., 11., 6., 3.]);
    
    let (sum_p, diam, leach, time) = (600.0, 2.0, 0.0, 1.0);

    let a = yasso20::get_a(&theta, &avg_t, sum_p, diam, leach);

    println!("--- Start Simulation ---");
    for year in 0..10 {
        init = yasso20::get_next_timestep(a, init, infall, time);
        println!("Year {year}: {:.3?}", init.as_slice());
    }

    let spin = yasso20::get_spin(a, infall);
    println!("Spin-up: {:.3?}", spin);

    // Benchmark
    println!("\n--- Start Benchmark (1M Iterationen) ---");
    let n = 1_000_000;
    let start = Instant::now();
    let mut result = SVector::<f64, 5>::zeros();
    
    for i in 1..=n {
        let factor = i as f64 / n as f64;
        // In Rust: .add_scalar(x) entspricht Julia's .+ x
        let avg_t2 = avg_t.add_scalar(factor);
        let a_dynamic = yasso20::get_a(&theta, &avg_t2, sum_p, diam, leach);
        result = yasso20::get_next_timestep(a_dynamic, result, infall, time);
    }
    
    let duration = start.elapsed();
    println!("Resultat: {:.5?}", result.as_slice());
    println!("Dauer: {:.4} Sekunden", duration.as_secs_f64());
}
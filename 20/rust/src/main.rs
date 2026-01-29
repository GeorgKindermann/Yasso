mod yasso;
use yasso::Yasso20;

fn main() {
    let mut model = Yasso20::new();

    let theta_vals: [f64; 35] = [
        0.51, 5.19, 0.13, 0.1, 0.5, 0., 1., 1., 0.99, 0., 0., 0., 0., 0., 0.163, 0., -0., 0., 0., 0.,
        0., 0.158, -0.002, 0.17, -0.005, 0.067, -0., -1.44, -2.0, -6.9, 0.0042, 0.0015, -2.55,
        1.24, 0.25,
    ];

    let mut init: [f64; 5] = [0.0; 5];
    let infall: [f64; 5] = [0.5, 0.1, 0.1, 0.2, 0.0];
    let time = 1.0;
    let avg_t: [f64; 12] = [
        -2.0, 2.0, 5.0, 10.0, 15.0, 20.0, 22.0, 21.0, 16.0, 11.0, 6.0, 3.0,
    ];
    let sum_p = 600.0;
    let diam = 2.0;
    let leach = 0.0;

    model.set_theta(theta_vals);
    model.set_timespan(time);
    
    model.set_clim_size_leach(&avg_t, sum_p, diam, leach);

    println!("Simulation (first 10 Years):");
    for year in 0..10 {
        let mut next_res = [0.0; 5];
        model.get_next_timestep(&init, &infall, &mut next_res);
        init = next_res;
        println!("{} {:?}", year, init);
    }

    let mut spin_res = [0.0; 5];
    model.get_spin(&infall, &mut spin_res, 700);
    println!("*: {:?}", spin_res);

    let n = 1_000_000;
    let mut result = [0.0; 5];
    let start = std::time::Instant::now();

    for i in 1..=n {
        let factor = i as f64 / n as f64;
        
        let mut current_infall = [0.0; 5];
        for k in 0..5 { current_infall[k] = infall[k] * factor; }
        
        let mut current_avg_t = [0.0; 12];
        for k in 0..12 { current_avg_t[k] = avg_t[k] + factor; }

        model.set_clim_size_leach(&current_avg_t, sum_p, diam, leach);
        
        let mut next_step = [0.0; 5];
        model.get_next_timestep(&result, &current_infall, &mut next_step);
        result = next_step;
    }

    let duration = start.elapsed();
    println!("Time: {:.4} Sekunden", duration.as_secs_f64());
    println!("Result: {:?}", result);
}

mod yasso20;
use yasso20::Yasso20;
use std::time::Instant;

fn main() {
    let theta = [0.51,5.19,0.13,0.1,0.5,0.,1.,1.,0.99,0.,0.,0.,0.,0.,0.163,0.,-0.,0.,0.,0.,0.,0.158,-0.002,0.17,-0.005,0.067,-0.,-1.44,-2.0,-6.9,0.0042,0.0015,-2.55,1.24,0.25];
    let avg_t = [-2., 2., 5., 10., 15., 20., 22., 21., 16., 11., 6., 3.];
    let mut init = [0.0; 5];
    let infall = [0.5, 0.1, 0.1, 0.2, 0.0];
    
    let mut yasso = Yasso20::new(theta);
    yasso.set_clim_size_leach(&avg_t, 600.0, 2.0, 0.0);
    for year in 0..10 {
      let mut next_init = [0.0; 5];
      yasso.get_next_timestep(&init, &infall, &mut next_init);
      init = next_init;
      println!("{} [{:.3}, {:.3}, {:.3}, {:.3}, {:.3}]", 
              year, init[0], init[1], init[2], init[3], init[4]);
    }

    let spin = yasso.get_spin(&infall);
    println!("Spin-up: {:.3?}", spin);
    
    let n = 1_000_000;
    let start = Instant::now();
  for i in 1..=n {
    let mut infall2 = [0.0; 5];
    let mut avg_t2 = [0.0; 12];
    let factor = (i as f64) / (n as f64);
    for j in 0..5 { infall2[j] = infall[j] * factor; }
    for j in 0..12 { avg_t2[j] = avg_t[j] + factor; }
    
    yasso.set_clim_size_leach(&avg_t2, 600.0, 2.0, 0.0);
    
    // temporären Speicher für das Ergebnis
    let mut next_init = [0.0; 5];
    yasso.get_next_timestep(&init, &infall2, &mut next_init);
    init = next_init;
}
    println!("Time: {:.6}s", start.elapsed().as_secs_f64());
    println!("Result: {:.5?}", init);
}
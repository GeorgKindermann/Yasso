use nalgebra::{SMatrix, SVector};

#[inline(always)]
pub fn exp5(a: SMatrix<f64, 5, 5>, q: usize) -> SMatrix<f64, 5, 5> {
    // Norm-Berechnung
    let s: f64 = a.iter().map(|x| x * x).sum();
    let p = s.sqrt();
    
    // Scaling-Phase
    let mut normiter = 2.0;
    let mut j = 2;
    while p > normiter {
        normiter *= 2.0;
        j += 1;
    }

    let c = a / normiter;
    let mut res = SMatrix::<f64, 5, 5>::identity();
    res += c;
    
    // Taylor-Reihe für q Terme
    let mut d = c;
    for i in 2..q {
        d *= c;
        // Fakultät iterativ anwenden: d = d / i
        d /= i as f64;
        res += d;
    }

    // Squaring-Phase
    for _ in 1..j {
        res = res * res;
    }
    
    res
}

#[allow(dead_code)]
#[inline(always)]
pub fn exp5_fix(a: SMatrix<f64, 5, 5>) -> SMatrix<f64, 5, 5> {
    let s: f64 = a.iter().map(|x| x * x).sum();
    let p = s.sqrt();
    
    let mut normiter = 2.0;
    let mut j = 2;
    while p > normiter {
        normiter *= 2.0;
        j += 1;
    }

    let c = a / normiter;
    let mut res = SMatrix::<f64, 5, 5>::identity() + c;
    
    let mut d = c * c;
    res += d * 0.5;
    d *= c; res += d * 0.16666666666666666;
    d *= c; res += d * 0.041666666666666664;
    d *= c; res += d * 0.008333333333333333;
    d *= c; res += d * 0.001388888888888889;
    d *= c; res += d * 0.0001984126984126984;

    for _ in 2..=j {
        res = res * res;
    }
    res
}

#[inline(always)]
pub fn get_a(theta: &SVector<f64, 35>, avg_t: &SVector<f64, 12>, sum_p: f64, diam: f64, leach: f64) -> SMatrix<f64, 5, 5> {
    let m3 = sum_p / 1000.0;

    let mut tem_sum = 0.0;
    let mut tem_n_sum = 0.0;
    let mut tem_h_sum = 0.0;

    for &t in avg_t.iter() {
        let t2 = t * t;
        tem_sum += (theta[21] * t + theta[22] * t2).exp();
        tem_n_sum += (theta[23] * t + theta[24] * t2).exp();
        tem_h_sum += (theta[25] * t + theta[26] * t2).exp();
    }

    let tem = (tem_sum / 12.0) * (1.0 - (theta[27] * m3).exp());
    let tem_n = (tem_n_sum / 12.0) * (1.0 - (theta[28] * m3).exp());
    let tem_h = (tem_h_sum / 12.0) * (1.0 - (theta[29] * m3).exp());

    let size_dep = if diam > 0.0 {
        (1.0 + theta[32] * diam + theta[33] * (diam * diam)).powf(theta[34]).min(1.0)
    } else { 1.0 };

    let mut a = SMatrix::<f64, 5, 5>::zeros();
    let k = tem * size_dep;
    
    a[(0, 0)] = theta[0] * k;
    a[(1, 1)] = theta[1] * k;
    a[(2, 2)] = theta[2] * k;
    a[(3, 3)] = theta[3] * tem_n * size_dep;
    a[(4, 4)] = theta[31] * tem_h;

    let d_abs = [a[(0,0)].abs(), a[(1,1)].abs(), a[(2,2)].abs(), a[(3,3)].abs()];

    let mut idx = 4; // Entspricht Julia Index 5
    for i in 0..4 {
        for j in 0..4 {
            if i != j {
                a[(i, j)] = theta[idx] * d_abs[j];
                idx += 1;
            }
        }
    }

    for j in 0..4 { a[(4, j)] = theta[30] * d_abs[j]; }

    if leach < 0.0 {
        let lm3 = leach * m3;
        for j in 0..4 { a[(j, j)] += lm3; }
    }
    a
}

#[inline(always)]
pub fn get_next_timestep(a: SMatrix<f64, 5, 5>, init: SVector<f64, 5>, infall: SVector<f64, 5>, t: f64) -> SVector<f64, 5> {
    let exp_at = exp5(a * t, 6);
    let rhs = exp_at * (a * init + infall) - infall;
    a.lu().solve(&rhs).expect("Matrix inversion failed")
}

pub fn get_theta(mut theta_vec: Vec<f64>) -> SVector<f64, 35> {
    // Julia 1,2,3,4, 32, 35 -> Rust 0, 1, 2, 3, 31, 34
    let tabs = [0, 1, 2, 3, 31, 34]; 
    for &i in &tabs {
        theta_vec[i] = -theta_vec[i].abs();
    }
    SVector::<f64, 35>::from_vec(theta_vec)
}

#[inline(always)]
pub fn get_spin(a: SMatrix<f64, 5, 5>, infall: SVector<f64, 5>) -> SVector<f64, 5> {
    a.lu().solve(&infall).map(|res| res * -1.0).expect("Matrix inversion failed")
}
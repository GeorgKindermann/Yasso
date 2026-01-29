pub struct Yasso20 {
    theta: [f64; 35],
    a: [f64; 25],
    a_decomp: [f64; 25],
    mexp_at: [f64; 25],
    timespan: f64,
    taylor_terms: usize,
    tol: f64,
    a_need_to_be_decomp: bool,
    a_need_to_be_expo: bool,
    no_decomposition: bool,
}

impl Yasso20 {
    pub fn new() -> Self {
        Self {
            theta: [0.0; 35],
            a: [0.0; 25],
            a_decomp: [0.0; 25],
            mexp_at: [0.0; 25],
            timespan: 1.0,
            taylor_terms: 11,
            tol: 1e-12,
            a_need_to_be_decomp: true,
            a_need_to_be_expo: true,
            no_decomposition: true,
        }
    }

    pub fn set_theta(&mut self, theta_input: [f64; 35]) {
        self.theta = theta_input;
        self.theta[31] = -self.theta[31].abs();
        self.theta[34] = -self.theta[34].abs();
        for i in 0..4 {
            self.theta[i] = -self.theta[i].abs();
        }
    }

    pub fn set_timespan(&mut self, timespan: f64) {
        self.timespan = timespan;
        self.a_need_to_be_expo = true;
    }

    pub fn set_clim_size_leach(&mut self, avg_t: &[f64; 12], sum_p: f64, diam: f64, leach: f64) {
        let mut tem = 0.0;
        let mut tem_n = 0.0;
        let mut tem_h = 0.0;
        let m3 = sum_p / 1000.0;

        for &t in avg_t {
            let t2 = t * t;
            tem += (self.theta[21] * t + self.theta[22] * t2).exp();
            tem_n += (self.theta[23] * t + self.theta[24] * t2).exp();
            tem_h += (self.theta[25] * t + self.theta[26] * t2).exp();
        }
        tem /= 12.0; tem_n /= 12.0; tem_h /= 12.0;

        tem *= 1.0 - (self.theta[27] * m3).exp();
        tem_n *= 1.0 - (self.theta[28] * m3).exp();
        tem_h *= 1.0 - (self.theta[29] * m3).exp();

        if tem < self.tol {
            self.no_decomposition = true;
            return;
        }
        self.no_decomposition = false;

        self.a.fill(0.0);
        if diam > 0.0 {
            let size_dep = (1.0 + self.theta[32] * diam + self.theta[33] * diam * diam)
                .powf(self.theta[34]).min(1.0);
            for i in 0..3 { self.a[i * 6] = self.theta[i] * tem * size_dep; }
            self.a[3 * 6] = self.theta[3] * tem_n * size_dep;
        } else {
            for i in 0..3 { self.a[i * 6] = self.theta[i] * tem; }
            self.a[3 * 6] = self.theta[3] * tem_n;
        }

        self.a[24] = self.theta[31] * tem_h;
        let d_abs = [self.a[0].abs(), self.a[6].abs(), self.a[12].abs(), self.a[18].abs()];

        let mut idx = 4;
        for i in 0..4 {
            for j in 0..4 {
                if i != j {
                    self.a[i * 5 + j] = self.theta[idx] * d_abs[j];
                    idx += 1;
                }
            }
        }
        // Mass flows AWEN -> H
        for i in 0..4 { self.a[20 + i] = self.theta[30] * d_abs[i]; }

        if leach < 0.0 {
            let aux = leach * m3;
            for i in 0..4 { self.a[i * 6] += aux; }
        }

        self.a_need_to_be_decomp = true;
        self.a_need_to_be_expo = true;
    }

    pub fn get_spin(&mut self, infall: &[f64; 5], result: &mut [f64; 5], years: i32) {
        if self.no_decomposition {
            for i in 0..5 { result[i] = years as f64 * infall[i]; }
        } else {
            if self.a_need_to_be_decomp {
                crout::<5>(&self.a, &mut self.a_decomp);
                self.a_need_to_be_decomp = false;
            }
            solve_crout::<5>(&self.a_decomp, infall, result);
            for i in 0..5 { result[i] *= -1.0; }
        }
    }

    pub fn get_next_timestep(&mut self, init: &[f64; 5], infall: &[f64; 5], result: &mut [f64; 5]) {
        if self.no_decomposition {
            for i in 0..5 { result[i] += infall[i]; }
            return;
        }

        let mut z1 = [0.0; 5];
        matmul_vec_5(&self.a, init, &mut z1);
        for i in 0..5 { z1[i] += infall[i]; }

        if self.a_need_to_be_expo {
            let mut at = [0.0; 25];
            for i in 0..25 { at[i] = self.a[i] * self.timespan; }
            self.mexp_at = matrixexp_5(&at, self.taylor_terms);
            self.a_need_to_be_expo = false;
        }

        let mut z2 = [0.0; 5];
        matmul_vec_5(&self.mexp_at, &z1, &mut z2);
        for i in 0..5 { z2[i] -= infall[i]; }

        if self.a_need_to_be_decomp {
            crout::<5>(&self.a, &mut self.a_decomp);
            self.a_need_to_be_decomp = false;
        }
        solve_crout::<5>(&self.a_decomp, &z2, result);
    }
}

fn crout<const D: usize>(a: &[f64; 25], decomp: &mut [f64; 25]) {
    for k in 0..D {
        for i in k..D {
            let mut sum = 0.0;
            for p in 0..k { sum += decomp[i * D + p] * decomp[p * D + k]; }
            decomp[i * D + k] = a[i * D + k] - sum;
        }
        for j in k + 1..D {
            let mut sum = 0.0;
            for p in 0..k { sum += decomp[k * D + p] * decomp[p * D + j]; }
            decomp[k * D + j] = (a[k * D + j] - sum) / decomp[k * D + k];
        }
    }
}

fn solve_crout<const D: usize>(decomp: &[f64; 25], b: &[f64; 5], x: &mut [f64; 5]) {
    let mut y = [0.0; 5];
    for i in 0..D {
        let mut sum = 0.0;
        for k in 0..i { sum += decomp[i * D + k] * y[k]; }
        y[i] = (b[i] - sum) / decomp[i * D + i];
    }
    for i in (0..D).rev() {
        let mut sum = 0.0;
        for k in i + 1..D { sum += decomp[i * D + k] * x[k]; }
        x[i] = y[i] - sum;
    }
}

fn matmul_5(a: &[f64; 25], b: &[f64; 25], res: &mut [f64; 25]) {
    let mut tmp = [0.0; 25];
    for m in 0..5 {
        for p in 0..5 {
            let mut sum = 0.0;
            for n in 0..5 { sum += a[m * 5 + n] * b[n * 5 + p]; }
            tmp[m * 5 + p] = sum;
        }
    }
    res.copy_from_slice(&tmp);
}

fn matmul_vec_5(a: &[f64; 25], b: &[f64; 5], res: &mut [f64; 5]) {
    for m in 0..5 {
        let mut sum = 0.0;
        for n in 0..5 { sum += a[m * 5 + n] * b[n]; }
        res[m] = sum;
    }
}

fn matrixexp_5(a: &[f64; 25], q: usize) -> [f64; 25] {
    let mut b = [0.0; 25];
    for i in 0..5 { b[i * 6] = 1.0; }
    let norm = a.iter().map(|x| x * x).sum::<f64>().sqrt();
    let mut normiter = 2.0;
    let mut j = 2;
    while norm > normiter { normiter *= 2.0; j += 1; }
    let mut c = [0.0; 25];
    for i in 0..25 { c[i] = a[i] / normiter; }
    let mut d = c;
    for i in 0..25 { b[i] += c[i]; }
    for i in 2..q {
        let mut next_d = [0.0; 25];
        matmul_5(&c, &d, &mut next_d);
        d = next_d;
        for k in 0..25 { d[k] /= i as f64; b[k] += d[k]; }
    }
    for _ in 1..j {
        let mut next_b = [0.0; 25];
        matmul_5(&b, &b, &mut next_b);
        b = next_b;
    }
    b
}

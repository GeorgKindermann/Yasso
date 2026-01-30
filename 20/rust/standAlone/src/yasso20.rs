pub struct Yasso20 {
    theta: [f64; 35],
    a: [f64; 25],
    a_decomp: [f64; 25],
    mexp_at: [f64; 25],
    timespan: f64,
    taylor_terms: usize,
    no_decomposition: bool,
    a_need_to_decomp: bool,
    a_need_to_expo: bool,
}

impl Yasso20 {
    pub fn new(atheta: [f64; 35]) -> Self {
        let mut y = Yasso20 {
            theta: atheta,
            a: [0.0; 25],
            a_decomp: [0.0; 25],
            mexp_at: [0.0; 25],
            timespan: 1.0,
            taylor_terms: 7,
            no_decomposition: true,
            a_need_to_decomp: true,
            a_need_to_expo: true,
        };
        y.set_theta(atheta);
        y
    }

    

pub fn set_theta(&mut self, mut atheta: [f64; 35]) {
        atheta[31] = -atheta[31].abs();
        atheta[34] = -atheta[34].abs();
        for i in 0..4 { atheta[i] = -atheta[i].abs(); }
        self.theta = atheta;
    }

    #[inline(always)]
    pub fn set_clim_size_leach(&mut self, avg_t: &[f64; 12], sum_p: f64, diam: f64, leach: f64) {
        let mut tem = 0.0;
        let mut tem_n = 0.0;
        let mut tem_h = 0.0;
        let m3 = sum_p / 1000.0;

        for i in 0..12 {
            let t = avg_t[i];
            let t2 = t * t;
            tem += (self.theta[21] * t + self.theta[22] * t2).exp();
            tem_n += (self.theta[23] * t + self.theta[24] * t2).exp();
            tem_h += (self.theta[25] * t + self.theta[26] * t2).exp();
        }
        tem /= 12.0; tem_n /= 12.0; tem_h /= 12.0;
        tem *= 1.0 - (self.theta[27] * m3).exp();
        tem_n *= 1.0 - (self.theta[28] * m3).exp();
        tem_h *= 1.0 - (self.theta[29] * m3).exp();

        if tem < 1e-12 {
            self.no_decomposition = true;
            return;
        }
        self.no_decomposition = false;

        self.a = [0.0; 25];
        let size_dep = if diam > 0.0 {
            (1.0 + self.theta[32] * diam + self.theta[33] * diam * diam)
                .powf(self.theta[34]).min(1.0)
        } else { 1.0 };

        for i in 0..3 { self.a[i * 6] = self.theta[i] * tem * size_dep; }
        self.a[3 * 6] = self.theta[3] * tem_n * size_dep;
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
        for i in 0..4 { self.a[20 + i] = self.theta[30] * d_abs[i]; }
        if leach < 0.0 {
            let aux = leach * m3;
            for i in 0..4 { self.a[6 * i] += aux; }
        }
        self.a_need_to_decomp = true;
        self.a_need_to_expo = true;
    }

    #[inline(always)]
    pub fn get_next_timestep(&mut self, init: &[f64; 5], infall: &[f64; 5], result: &mut [f64; 5]) {
        if self.no_decomposition {
            for i in 0..5 { result[i] += infall[i]; }
            return;
        }

        let mut z1 = [0.0; 5];
        matmul::<5, 5, 1>(&self.a, init, &mut z1);
        for i in 0..5 { z1[i] += infall[i]; }

        if self.a_need_to_expo {
            let mut at = [0.0; 25];
            for i in 0..25 { at[i] = self.a[i] * self.timespan; }
            self.mexp_at = [0.0; 25];
            matrixexp::<5>(&at, &mut self.mexp_at, self.taylor_terms);
            self.a_need_to_expo = false;
        }

        let mut z2 = [0.0; 5];
        matmul::<5, 5, 1>(&self.mexp_at, &z1, &mut z2);
        for i in 0..5 { z2[i] -= infall[i]; }

        if self.a_need_to_decomp {
            crout::<5>(&self.a, &mut self.a_decomp);
            self.a_need_to_decomp = false;
        }
        solve_crout::<5>(&self.a_decomp, &z2, result);
    }

    pub fn get_spin(&mut self, infall: &[f64; 5]) -> [f64; 5] {
        let mut result = [0.0; 5];
        if self.a_need_to_decomp {
            crout::<5>(&self.a, &mut self.a_decomp);
            self.a_need_to_decomp = false;
        }
        solve_crout::<5>(&self.a_decomp, infall, &mut result);
        for i in 0..5 { result[i] *= -1.0; }
        result
    }
}

// --- Hilfsfunktionen ---

#[inline(always)]
fn matmul<const M: usize, const N: usize, const P: usize>(a: &[f64], b: &[f64], c: &mut [f64]) {
    for p in 0..P {
        for m in 0..M {
            let mut sum = a[m * N] * b[p];
            for n in 1..N {
                sum += a[m * N + n] * b[n * P + p];
            }
            c[m * P + p] = sum;
        }
    }
}

#[inline(always)]
fn crout<const D: usize>(a: &[f64; 25], decomp: &mut [f64; 25]) {
    for k in 0..D {
        for i in k..D {
            let mut sum = 0.0;
            for p in 0..k { sum += decomp[i * D + p] * decomp[p * D + k]; }
            decomp[i * D + k] = a[i * D + k] - sum;
        }
        for j in (k + 1)..D {
            let mut sum = 0.0;
            for p in 0..k { sum += decomp[k * D + p] * decomp[p * D + j]; }
            decomp[k * D + j] = (a[k * D + j] - sum) / decomp[k * D + k];
        }
    }
}

#[inline(always)]
fn solve_crout<const D: usize>(decomp: &[f64; 25], b: &[f64; 5], x: &mut [f64; 5]) {
    let mut y = [0.0; 5];
    for i in 0..D {
        let mut sum = 0.0;
        for k in 0..i { sum += decomp[i * D + k] * y[k]; }
        y[i] = (b[i] - sum) / decomp[i * D + i];
    }
    for i in (0..D).rev() {
        let mut sum = 0.0;
        for k in (i + 1)..D { sum += decomp[i * D + k] * x[k]; }
        x[i] = y[i] - sum;
    }
}

#[inline(always)]
fn matrixexp<const M: usize>(a: &[f64; 25], b: &mut [f64; 25], q: usize) {
    for i in 0..M { b[i * (M + 1)] = 1.0; }
    let mut normiter = 2.0;
    let mut j = 2;
    let mut p_sq = 0.0;
    for x in a { p_sq += x * x; }
    let p = p_sq.sqrt();
    while p > normiter {
        normiter *= 2.0;
        j += 1;
    }
    let mut c = [0.0; 25];
    for i in 0..25 { c[i] = a[i] / normiter; }
    for i in 0..25 { b[i] += c[i]; }
    let mut d = c;
    for i in 2..q {
        let mut next_d = [0.0; 25];
        matmul::<M, M, M>(&c, &d, &mut next_d);
        d = next_d;
        let inv_i = 1.0 / (i as f64);
        for j in 0..25 { d[j] *= inv_i; }
        for j in 0..25 { b[j] += d[j]; }
    }
    for _ in 1..j {
        let mut next_b = [0.0; 25];
        matmul::<M, M, M>(b, b, &mut next_b);
        *b = next_b;
    }
}
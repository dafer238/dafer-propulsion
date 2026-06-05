use std::f64::consts::PI;

pub trait GasModel {
    fn gamma(&self) -> f64;

    fn prandtl_meyer(&self, m: f64) -> f64;
    fn inverse_prandtl_meyer(&self, nu: f64) -> f64;
}

pub struct Air {
    gamma: f64,
}

impl Air {
    pub fn new(gamma: f64) -> Self {
        Self { gamma }
    }
}

impl GasModel for Air {
    fn gamma(&self) -> f64 {
        self.gamma
    }

    fn prandtl_meyer(&self, m: f64) -> f64 {
        let g = self.gamma;
        let a = (g + 1.0) / (g - 1.0);
        (a.sqrt() * ( ((g - 1.0)/(g + 1.0) * (m*m - 1.0)).sqrt() ).atan()
        - (m*m - 1.0).sqrt().atan()) * 180.0 / PI
    }

    fn inverse_prandtl_meyer(&self, nu: f64) -> f64 {
        // very rough bisection solve
        let mut m = 2.0;
        for _ in 0..50 {
            let nu_m = self.prandtl_meyer(m);
            let err = nu_m - nu;
            m -= err * 0.01;
            if m < 1.0 { m = 1.01; }
        }
        m
    }
}

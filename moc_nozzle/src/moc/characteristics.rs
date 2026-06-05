use crate::core::state::FlowState;

pub struct Invariants {
    pub k_plus: f64,
    pub k_minus: f64,
}

pub fn invariants(s: FlowState) -> Invariants {
    Invariants {
        k_plus: s.theta + s.nu,
        k_minus: s.theta - s.nu,
    }
}

pub fn from_invariants(kp: f64, km: f64) -> FlowState {
    FlowState {
        theta: (kp + km) / 2.0,
        nu: (kp - km) / 2.0,
        m: 2.0, // placeholder, overwritten later
    }
}

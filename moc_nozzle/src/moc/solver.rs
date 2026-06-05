use crate::moc::node::Node;
use crate::core::state::FlowState;
use crate::moc::characteristics::{invariants, from_invariants};

pub struct SimpleMocSolver {
    pub n: usize,
    pub gamma: f64,
    pub nodes: Vec<Node>,
}

impl SimpleMocSolver {
    pub fn new(n: usize, gamma: f64) -> Self {
        Self {
            n,
            gamma,
            nodes: Vec::new(),
        }
    }

    pub fn nodes(&self) -> &[Node] {
        &self.nodes
    }

    pub fn initialize(&mut self) {
        for i in 0..self.n {
            let theta = i as f64 / self.n as f64 * 0.5;
            let nu = theta;

            self.nodes.push(Node {
                x: i as f64,
                y: 0.0,
                state: FlowState { m: 2.0, theta, nu },
            });
        }
    }

    pub fn step(&mut self) {
        if self.nodes.is_empty() {
            self.initialize();
        }

        let mut new_nodes = Vec::new();

        for i in 1..self.nodes.len() {
            let a = &self.nodes[i - 1];
            let b = &self.nodes[i];

            let ka = invariants(a.state);
            let kb = invariants(b.state);

            let state = from_invariants(ka.k_plus, kb.k_minus);

            new_nodes.push(Node {
                x: (a.x + b.x) / 2.0,
                y: (a.y + b.y) / 2.0 + 0.1,
                state,
            });
        }

        self.nodes = new_nodes;
    }
}

use crate::core::gas::GasModel;
use crate::{moc::solver::SimpleMocSolver, solver::config::NozzleConfig};

pub struct NozzleSolver<G: GasModel> {
    pub gas: G,
    pub solver: SimpleMocSolver,
    pub config: NozzleConfig,
}

impl<G: GasModel> NozzleSolver<G> {
    pub fn new(gas: G, solver: SimpleMocSolver, config: NozzleConfig) -> Self {
        Self {
            gas,
            solver,
            config,
        }
    }

    pub fn run(&mut self) {
        self.solver.initialize();

        for _ in 0..5 {
            self.solver.step();
        }
    }
}

mod core;
mod geometry;
mod moc;
mod solver;
mod utils;

use core::gas::Air;
use moc::solver::SimpleMocSolver;
use solver::config::NozzleConfig;
use solver::nozzle::NozzleSolver;

fn main() {
    let gas = Air::new(1.4);

    let config = NozzleConfig {
        gamma: 1.4,
        ae_at: 10.0,
        n_points: 50,
    };

    let mut solver = NozzleSolver::new(
        gas,
        SimpleMocSolver::new(config.n_points, config.gamma),
        config,
    );

    solver.run();

    let nodes = solver.solver.nodes();

    println!("Generated {} nodes", nodes.len());

    for n in nodes.iter().take(50) {
        println!("{:?}", n);
    }
}

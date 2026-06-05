use crate::core::state::FlowState;

#[derive(Clone, Debug)]
pub struct Node {
    pub x: f64,
    pub y: f64,
    pub state: FlowState,
}

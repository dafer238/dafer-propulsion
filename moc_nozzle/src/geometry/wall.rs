use crate::moc::node::Node;

pub struct NozzleWall {
    pub points: Vec<(f64, f64)>,
}

pub fn extract_wall(nodes: &[Node]) -> NozzleWall {
    let mut pts = Vec::new();

    for n in nodes {
        if n.y >= 0.5 {
            pts.push((n.x, n.y));
        }
    }

    NozzleWall { points: pts }
}

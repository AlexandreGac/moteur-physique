mod rigid_body;
mod solver;
mod collision;
mod constraints;

use nalgebra::Point2;
use parry2d::shape::{ConvexPolygon, Shape};
use crate::rigid_body::RigidBody;
use crate::solver::SolverBuilder;

fn main() {
    let points1 = vec![
        Point2::new(0.0, 0.0),
        Point2::new(1.0, 0.0),
        Point2::new(1.0, 1.0),
        Point2::new(0.0, 1.0)
    ];
    let polygon1 = ConvexPolygon::from_convex_polyline(points1).unwrap();
    let rb1 = RigidBody::new(polygon1, 1.0);

    let points2 = vec![
        Point2::new(0.5, 0.9),
        Point2::new(1.0, 1.9),
        Point2::new(0.0, 1.9)
    ];
    let polygon2 = ConvexPolygon::from_convex_polyline(points2).unwrap();
    let rb2 = RigidBody::new(polygon2, 1.0);

    let points3 = vec![
        Point2::new(0.8, 0.8),
        Point2::new(1.8, 0.8),
        Point2::new(1.8, 1.8),
        Point2::new(0.8, 1.8)
    ];
    let polygon3 = ConvexPolygon::from_convex_polyline(points3).unwrap();
    let rb3 = RigidBody::new(polygon3, 1.0);

    let mut solver = SolverBuilder::new()
        .add_rigid_body(rb1)
        .add_rigid_body(rb2)
        .add_rigid_body(rb3)
        .build();

    solver.simulate();
}

mod rigid_body;
mod solver;
mod collision;
mod constraints;
mod utils;

use nalgebra::{Point2, Vector2};
use parry2d_f64::math::Real;
use parry2d_f64::shape::{ConvexPolygon, Shape};
use crate::rigid_body::RigidBody;
use crate::solver::{Solver, SolverBuilder};
use piston_window::*;
use piston_window::math::Scalar;
use crate::constraints::Constraint;

fn main() {
    let points_0 = vec![
        Point2::new(-5.5, 4.5),
        Point2::new(-4.5, 4.5),
        Point2::new(-4.5, 5.5),
        Point2::new(-5.5, 5.5)
    ];
    let polygon_0 = ConvexPolygon::from_convex_polyline(points_0).unwrap();
    let rb_0 = RigidBody::new(polygon_0, 1.0);

    let points_1 = vec![
        Point2::new(-5.5, 0.5),
        Point2::new(5.5, 0.5),
        Point2::new(5.5, 1.5),
        Point2::new(-5.5, 2.5)
    ];
    let polygon_1 = ConvexPolygon::from_convex_polyline(points_1).unwrap();
    let rb_1 = RigidBody::new(polygon_1, Real::INFINITY);

    let solver = SolverBuilder::new()
        .add_rigid_body(rb_0)
        .add_rigid_body(rb_1)
        .build();

    simulation_loop(solver);
}

fn simulation_loop(mut solver: Solver) {
    let (width, height) = (1200, 700);
    let mut window: PistonWindow = WindowSettings::new("Moteur physique", [width, height])
        .vsync(true)
        .exit_on_esc(true)
        .build()
        .unwrap();

    while let Some(event) = window.next() {
        solver.simulate();
        window.draw_2d(&event, |context, graphics, _device| {
            clear([0.4; 4], graphics);
            solver.display(context, graphics, width as Scalar, height as Scalar);
        });
    }
}

mod rigid_body;
mod solver;
mod collision;
mod constraints;
mod utils;

use nalgebra::{Isometry2, Point2, Vector2};
use parry2d_f64::math::Real;
use parry2d_f64::shape;
use parry2d_f64::shape::ConvexPolygon;
use crate::rigid_body::RigidBody;
use crate::solver::{Solver, SolverBuilder};
use piston_window::*;
use piston_window::math::Scalar;
use crate::constraints::Constraint;

fn main() {

    let shape_0 = shape::Ball::new(0.5);
    let mut rb_0 = RigidBody::new(shape_0, 1.0);
    rb_0.set_transform(Isometry2::translation(-5.0, 5.0));
    let pendulum_0 = Constraint::create_distance_constraint(
        0,
        Vector2::new(0.0, 0.0),
        Vector2::new(-1.0, 5.0)
    );

    let points_1 = vec![
        Point2::new(-1.5, 0.5),
        Point2::new(-0.5, 0.5),
        Point2::new(-0.5, 1.5),
        Point2::new(-1.0, 1.5)
    ];
    let polygon_1 = ConvexPolygon::from_convex_polyline(points_1).unwrap();
    let rb_1 = RigidBody::new(polygon_1, Real::INFINITY);

    let solver = SolverBuilder::new()
        .add_rigid_body(rb_0)
        .add_rigid_body(rb_1)
        .add_constraint(pendulum_0)
        .build();

/*
    let points_0 = vec![
        Point2::new(-5.5, 4.5),
        Point2::new(-4.5, 4.5),
        Point2::new(-4.5, 5.5),
        Point2::new(-5.5, 5.5)
    ];
    let polygon_0 = ConvexPolygon::from_convex_polyline(points_0).unwrap();
    let rb_0 = RigidBody::new(polygon_0, 1.0);
    let pendulum_0 = Constraint::create_distance_constraint(
        0,
        Vector2::new(0.0, 0.0),
        Vector2::new(-1.0, 5.0)
    );

    let points_1 = vec![
        Point2::new(-0.5, 0.5),
        Point2::new(0.5, 0.5),
        Point2::new(0.5, 1.5),
        Point2::new(-0.5, 1.5)
    ];
    let polygon_1 = ConvexPolygon::from_convex_polyline(points_1).unwrap();
    let rb_1 = RigidBody::new(polygon_1, 1.0);
    let pendulum_1 = Constraint::create_distance_constraint(
        1,
        Vector2::new(0.0, 0.0),
        Vector2::new(0.0, 5.0)
    );

    let points_2 = vec![
        Point2::new(0.5, 0.5),
        Point2::new(1.5, 0.5),
        Point2::new(1.5, 1.5),
        Point2::new(0.5, 1.5)
    ];
    let polygon_2 = ConvexPolygon::from_convex_polyline(points_2).unwrap();
    let rb_2 = RigidBody::new(polygon_2, 1.0);
    let pendulum_2 = Constraint::create_distance_constraint(
        2,
        Vector2::new(0.0, 0.0),
        Vector2::new(1.0, 5.0)
    );

    let solver = SolverBuilder::new()
        .add_rigid_body(rb_0)
        .add_rigid_body(rb_1)
        .add_rigid_body(rb_2)
        .add_constraint(pendulum_0)
        .add_constraint(pendulum_1)
        .add_constraint(pendulum_2)
        .build();
*/
/*
    let mut rb_0 = RigidBody::new(shape::Ball::new(0.5), 1.0);
    rb_0.set_transform(Isometry2::translation(0.0, 2.0));

    let points_2 = vec![
        Point2::new(-0.5001, 42.5),
        Point2::new(0.4999, 42.5),
        Point2::new(0.4999, 43.5),
        Point2::new(-0.5001, 43.5)
    ];
    let polygon_2 = ConvexPolygon::from_convex_polyline(points_2).unwrap();
    let rb_2 = RigidBody::new(polygon_2, 1.0);

    let points_1 = vec![
        Point2::new(-0.5, 0.5),
        Point2::new(0.5, 0.5),
        Point2::new(0.5, 1.5),
        Point2::new(-0.5, 1.5)
    ];
    let polygon_1 = ConvexPolygon::from_convex_polyline(points_1).unwrap();
    let rb_1 = RigidBody::new(polygon_1, Real::INFINITY);

    let solver = SolverBuilder::new()
        .add_rigid_body(rb_0)
        .add_rigid_body(rb_1)
        .add_rigid_body(rb_2)
        .build();
*/
/*
    let mut rb_2 = RigidBody::new(shape::Ball::new(0.5), 1.0);
    rb_2.set_transform(Isometry2::translation(-5.0, 10.0));

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
        .add_rigid_body(rb_2)
        .build();
*/
/*
    let mut rb_0 = RigidBody::new(shape::Ball::new(0.5), 1.0);
    rb_0.set_transform(Isometry2::translation(-5.0, 6.0));

    let points_1 = vec![
        Point2::new(-5.5, 0.5),
        Point2::new(-0.1, 0.5),
        Point2::new(-0.1, 1.5),
        Point2::new(-5.5, 3.5)
    ];
    let polygon_1 = ConvexPolygon::from_convex_polyline(points_1).unwrap();
    let rb_1 = RigidBody::new(polygon_1, Real::INFINITY);

    let points_2 = vec![
        Point2::new(0.1, 0.5),
        Point2::new(5.5, 0.5),
        Point2::new(5.5, 3.5),
        Point2::new(0.1, 1.5)
    ];
    let polygon_2 = ConvexPolygon::from_convex_polyline(points_2).unwrap();
    let rb_2 = RigidBody::new(polygon_2, Real::INFINITY);

    let points_3 = vec![
        Point2::new(4.5, 5.5),
        Point2::new(5.5, 5.5),
        Point2::new(5.5, 6.5),
        Point2::new(4.5, 6.5)
    ];
    let polygon_3 = ConvexPolygon::from_convex_polyline(points_3).unwrap();
    let rb_3 = RigidBody::new(polygon_3, 1.0);

    let points_4 = vec![
        Point2::new(4.5, 7.5),
        Point2::new(5.5, 7.5),
        Point2::new(5.5, 8.5),
        Point2::new(4.5, 8.5)
    ];
    let polygon_4 = ConvexPolygon::from_convex_polyline(points_4).unwrap();
    let rb_4 = RigidBody::new(polygon_4, 1.0);

    let points_5 = vec![
        Point2::new(-5.5, 7.5),
        Point2::new(-4.5, 7.5),
        Point2::new(-4.5, 8.5),
        Point2::new(-5.5, 8.5)
    ];
    let polygon_5 = ConvexPolygon::from_convex_polyline(points_5).unwrap();
    let rb_5 = RigidBody::new(polygon_5, 1.0);

    let points_6 = vec![
        Point2::new(-3.5, 16.5),
        Point2::new(-2.5, 16.5),
        Point2::new(-2.5, 17.5),
        Point2::new(-3.5, 17.5)
    ];
    let polygon_6 = ConvexPolygon::from_convex_polyline(points_6).unwrap();
    let rb_6 = RigidBody::new(polygon_6, 1.0);

    let points_7 = vec![
        Point2::new(2.5, 17.5),
        Point2::new(3.5, 17.5),
        Point2::new(3.5, 18.5),
        Point2::new(2.5, 18.5)
    ];
    let polygon_7 = ConvexPolygon::from_convex_polyline(points_7).unwrap();
    let rb_7 = RigidBody::new(polygon_7, 1.0);

    let points_8 = vec![
        Point2::new(-0.4, 2.0),
        Point2::new(0.6, 2.0),
        Point2::new(0.6, 3.0),
        Point2::new(-0.4, 3.0)
    ];
    let polygon_8 = ConvexPolygon::from_convex_polyline(points_8).unwrap();
    let rb_8 = RigidBody::new(polygon_8, 1.0);
    let pendulum_8 = Constraint::create_distance_constraint(
        8,
        Vector2::<Real>::new(0.5, 0.5),
        Vector2::<Real>::new(0.1, 5.0)
    );

    let points_9 = vec![
        Point2::new(-0.4, 22.0),
        Point2::new(0.6, 22.0),
        Point2::new(0.6, 23.0),
        Point2::new(-0.4, 23.0)
    ];
    let polygon_9 = ConvexPolygon::from_convex_polyline(points_9).unwrap();
    let rb_9 = RigidBody::new(polygon_9, 1.0);

    let solver = SolverBuilder::new()
        .add_rigid_body(rb_0)
        .add_rigid_body(rb_1)
        .add_rigid_body(rb_2)
        .add_rigid_body(rb_3)
        .add_rigid_body(rb_4)
        .add_rigid_body(rb_5)
        .add_rigid_body(rb_6)
        .add_rigid_body(rb_7)
        .add_rigid_body(rb_8)
        .add_rigid_body(rb_9)
        .add_constraint(pendulum_8)
        .build();
*/
    simulation_loop(solver);
}

fn simulation_loop(mut solver: Solver) {
    let (width, height) = (1200, 700);
    let mut window: PistonWindow = WindowSettings::new("Moteur physique", [width, height])
        .vsync(true)
        .exit_on_esc(true)
        .build()
        .unwrap();

    let mut debug = false;
    while let Some(event) = window.next() {
        if let Some(Button::Keyboard(Key::Space)) = event.press_args() {
            debug = !debug;
        }
        if !debug || true {
            debug = solver.simulate();
        }
        window.draw_2d(&event, |context, graphics, _device| {
            clear([0.4; 4], graphics);
            solver.display(context, graphics, width as Scalar, height as Scalar);
        });

    }
}

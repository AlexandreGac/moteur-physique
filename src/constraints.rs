use std::f64::consts::FRAC_PI_2;
use nalgebra::{Rotation2, RowDVector, Vector2};
use parry2d_f64::math::Real;
use piston_window::{Context, Graphics, line, rectangle};
use piston_window::math::Scalar;
use crate::rigid_body::{RigidBody, RigidBodyState};
use crate::solver::Index;

pub enum Constraint {
    Distance {
        index: Index,
        local: Vector2<Real>,
        pin: Vector2<Real>
    }
}

impl Constraint {

    pub fn create_distance_constraint(index: Index, local: Vector2<Real>, pin: Vector2<Real>) -> Constraint {
        Self::Distance {
            index, local, pin
        }
    }

    pub fn compute_jacobian(&self, states: &Vec<RigidBodyState>) -> RowDVector<Real> {
        let mut row = RowDVector::from(vec![0.0; 3 * states.len()]);
        match self {
            Constraint::Distance { index, local, pin } => {
                let state = &states[*index];
                let attachment: Vector2<Real> = Rotation2::new(state.position.z) * local;
                let position: Vector2<Real> = state.position.xy() + attachment;
                let orthogonal: Vector2<Real> = Rotation2::new(FRAC_PI_2) * attachment;

                row[3 * index] = 2.0 * (position.x - pin.x);
                row[3 * index + 1] = 2.0 * (position.y - pin.y);
                row[3 * index + 2] = 2.0 * orthogonal.dot(&(position - pin));
            }
            _ => todo!("Pas encore implémenté !")
        }

        row
    }

    pub fn compute_jacobian_derivative(&self, states: &Vec<RigidBodyState>) -> RowDVector<Real> {
        let mut row = RowDVector::from(vec![0.0; 3 * states.len()]);
        match self {
            Constraint::Distance { index, local, pin } => {
                let state = &states[*index];
                let attachment: Vector2<Real> = Rotation2::new(state.position.z) * local;
                let position: Vector2<Real> = state.position.xy() + attachment;
                let orthogonal: Vector2<Real> = Rotation2::new(FRAC_PI_2) * attachment;
                let velocity: Vector2<Real> = state.velocity.xy() + state.velocity.z * orthogonal;

                row[3 * index] = 2.0 * velocity.x;
                row[3 * index + 1] = 2.0 * velocity.y;
                row[3 * index + 2] = 2.0 * (orthogonal.dot(&velocity) - state.velocity.z * attachment.dot(&(position - pin)));
            }
            _ => todo!("Pas encore implémenté !")
        }

        row
    }

    pub fn display(&self, rigid_bodies: &Vec<RigidBody>, context: Context, graphics: &mut impl Graphics, width: Scalar, height: Scalar) {
        match self {
            Constraint::Distance { index, local, pin } => {
                let point_1 = Vector2::new(width as Real / 2.0 + pin.x * 100.0, height as Real - pin.y * 100.0);
                rectangle(
                    [1.0, 1.0, 1.0, 1.0],
                    [(point_1.x - 5.0) as Scalar, (point_1.y - 5.0) as Scalar, 10.0, 10.0],
                    context.transform,
                    graphics
                );

                let rotation = Rotation2::new(rigid_bodies[*index].position.z);
                let position = rotation * local + rigid_bodies[*index].position.xy();

                let point_2 = Vector2::new(width as Real / 2.0 + position.x * 100.0, height as Real - position.y * 100.0);
                rectangle(
                    [1.0, 1.0, 1.0, 1.0],
                    [(point_2.x - 5.0) as Scalar, (point_2.y - 5.0) as Scalar, 10.0, 10.0],
                    context.transform,
                    graphics
                );

                line([0.0, 0.0, 1.0, 1.0],
                     2.0,
                     [point_1.x as Scalar, point_1.y as Scalar, point_2.x as Scalar, point_2.y as Scalar],
                     context.transform,
                     graphics
                );
            },
            _ => {}
        }
    }
}
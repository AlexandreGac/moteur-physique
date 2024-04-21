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
    },
    Link {
        index_1: Index,
        index_2: Index,
        local_1: Vector2<Real>,
        local_2: Vector2<Real>
    }
}

impl Constraint {

    pub fn create_distance_constraint(index: Index, local: Vector2<Real>, pin: Vector2<Real>) -> Constraint {
        Self::Distance {
            index, local, pin
        }
    }

    pub fn create_link_constraint(index_1: Index, index_2: Index, local_1: Vector2<Real>, local_2: Vector2<Real>) -> Constraint {
        Self::Link {
            index_1, index_2, local_1, local_2
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
            Constraint::Link { index_1, index_2, local_1, local_2 } => {
                let state_1 = &states[*index_1];
                let state_2 = &states[*index_2];
                let attachment_1: Vector2<Real> = Rotation2::new(state_1.position.z) * local_1;
                let attachment_2: Vector2<Real> = Rotation2::new(state_2.position.z) * local_2;
                let position_1: Vector2<Real> = state_1.position.xy() + attachment_1;
                let position_2: Vector2<Real> = state_2.position.xy() + attachment_2;
                let direct_orthogonal = Rotation2::new(FRAC_PI_2);
                let orthogonal_1: Vector2<Real> = direct_orthogonal * attachment_1;
                let orthogonal_2: Vector2<Real> = direct_orthogonal * attachment_2;

                row[3 * index_1] = 2.0 * (position_1.x - position_2.x);
                row[3 * index_1 + 1] = 2.0 * (position_1.y - position_2.y);
                row[3 * index_1 + 2] = 2.0 * orthogonal_1.dot(&(position_1 - position_2));

                row[3 * index_2] = -2.0 * (position_1.x - position_2.x);
                row[3 * index_2 + 1] = -2.0 * (position_1.y - position_2.y);
                row[3 * index_2 + 2] = -2.0 * orthogonal_2.dot(&(position_1 - position_2));
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
            },
            Constraint::Link { index_1, index_2, local_1, local_2 } => {
                let state_1 = &states[*index_1];
                let state_2 = &states[*index_2];
                let attachment_1: Vector2<Real> = Rotation2::new(state_1.position.z) * local_1;
                let attachment_2: Vector2<Real> = Rotation2::new(state_2.position.z) * local_2;
                let position_1: Vector2<Real> = state_1.position.xy() + attachment_1;
                let position_2: Vector2<Real> = state_2.position.xy() + attachment_2;
                let direct_orthogonal = Rotation2::new(FRAC_PI_2);
                let orthogonal_1: Vector2<Real> = direct_orthogonal * attachment_1;
                let orthogonal_2: Vector2<Real> = direct_orthogonal * attachment_2;
                let velocity_1: Vector2<Real> = state_1.velocity.xy() + state_1.velocity.z * orthogonal_1;
                let velocity_2: Vector2<Real> = state_2.velocity.xy() + state_2.velocity.z * orthogonal_2;

                row[3 * index_1] = 2.0 * (velocity_1.x - velocity_2.x);
                row[3 * index_1 + 1] = 2.0 * (velocity_1.y - velocity_2.y);
                row[3 * index_1 + 2] = 2.0 * (orthogonal_1.dot(&(velocity_1 - velocity_2)) - state_1.velocity.z * attachment_1.dot(&(position_1 - position_2)));

                row[3 * index_2] = -2.0 * (velocity_1.x - velocity_2.x);
                row[3 * index_2 + 1] = -2.0 * (velocity_1.y - velocity_2.y);
                row[3 * index_2 + 2] = -2.0 * (orthogonal_2.dot(&(velocity_1 - velocity_2)) - state_2.velocity.z * attachment_2.dot(&(position_1 - position_2)));
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

                line([0.7, 0.7, 0.7, 1.0],
                     2.0,
                     [point_1.x as Scalar, point_1.y as Scalar, point_2.x as Scalar, point_2.y as Scalar],
                     context.transform,
                     graphics
                );
            },
            Constraint::Link { index_1, index_2, local_1, local_2 } => {
                let rotation_1 = Rotation2::new(rigid_bodies[*index_1].position.z);
                let position_1 = rotation_1 * local_1 + rigid_bodies[*index_1].position.xy();
                let point_1 = Vector2::new(width as Real / 2.0 + position_1.x * 100.0, height as Real - position_1.y * 100.0);
                rectangle(
                    [1.0, 1.0, 1.0, 1.0],
                    [(point_1.x - 5.0) as Scalar, (point_1.y - 5.0) as Scalar, 10.0, 10.0],
                    context.transform,
                    graphics
                );

                let rotation_2 = Rotation2::new(rigid_bodies[*index_2].position.z);
                let position_2 = rotation_2 * local_2 + rigid_bodies[*index_2].position.xy();
                let point_2 = Vector2::new(width as Real / 2.0 + position_2.x * 100.0, height as Real - position_2.y * 100.0);
                rectangle(
                    [1.0, 1.0, 1.0, 1.0],
                    [(point_2.x - 5.0) as Scalar, (point_2.y - 5.0) as Scalar, 10.0, 10.0],
                    context.transform,
                    graphics
                );

                line([0.7, 0.7, 0.7, 1.0],
                     2.0,
                     [point_1.x as Scalar, point_1.y as Scalar, point_2.x as Scalar, point_2.y as Scalar],
                     context.transform,
                     graphics
                );
            }
            _ => {}
        }
    }
}
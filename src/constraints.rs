use nalgebra::{RowDVector, Vector2};
use parry2d::math::Real;
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

    pub fn compute_jacobian(&self) -> RowDVector<Real> {
        RowDVector::<Real>::from(vec![1.0, 2.0, 3.0, 4.0, 5.0])
    }

    pub fn compute_jacobian_derivative(&self) -> RowDVector<Real> {
        RowDVector::<Real>::from(vec![1.0, 2.0, 3.0, 4.0, 5.0])
    }
}
use std::f64::consts::FRAC_PI_2;
use nalgebra::{RowDVector, Vector2, Rotation2};
use parry2d_f64::math::Real;
use parry2d_f64::query::{ContactManifold, DefaultQueryDispatcher, PersistentQueryDispatcher, TrackedContact};
use piston_window::math::Scalar;
use piston_window::{Context, Graphics, rectangle};
use crate::rigid_body::{RigidBody, RigidBodyState};
use crate::solver::Index;
use crate::utils::{COLLISION_EPS, CONTACT_VELOCITY_EPS};

#[derive(PartialEq, Clone, Debug)]
pub enum CollisionType {
    Penetration, Separation, Contact
}

#[derive(Clone, Debug)]
pub struct Collision {
    pub kind: CollisionType,
    pub index_1: Index,
    pub index_2: Index,
    pub point_1: Vector2<Real>,
    pub point_2: Vector2<Real>,
    pub normal_12: Vector2<Real>
}

impl Collision {

    pub fn update_collision_type(&mut self, rigid_bodies: &Vec<RigidBody>) {
        let rigid_body_1 = &rigid_bodies[self.index_1];
        let rigid_body_2 = &rigid_bodies[self.index_2];
        let direct_orthogonal = Rotation2::new(FRAC_PI_2);
        let velocity_1: Vector2<Real> = rigid_body_1.velocity.xy() + rigid_body_1.velocity.z * (direct_orthogonal * self.point_1);
        let velocity_2: Vector2<Real> = rigid_body_2.velocity.xy() + rigid_body_2.velocity.z * (direct_orthogonal * self.point_2);

        let relative_velocity_21 = (velocity_2 - velocity_1).dot(&self.normal_12);
        if relative_velocity_21.abs() < CONTACT_VELOCITY_EPS {
            self.kind = CollisionType::Contact;
        }
        else if relative_velocity_21 < 0.0 {
            self.kind = CollisionType::Penetration;
        }
        else {
            self.kind = CollisionType::Separation;
        }
    }

    pub fn get_restitution_coefficient(&self) -> Real {
        0.5
    }

    pub fn compute_jacobian(&self, rigid_bodies: &Vec<RigidBody>) -> RowDVector<Real> {
        let mut row = RowDVector::from(vec![0.0; 3 * rigid_bodies.len()]);
        let direct_orthogonal = Rotation2::new(FRAC_PI_2);
        row[3 * self.index_1] = -self.normal_12.x;
        row[3 * self.index_1 + 1] = -self.normal_12.y;
        row[3 * self.index_1 + 2] = -(direct_orthogonal * self.point_1).dot(&self.normal_12);

        row[3 * self.index_2] = self.normal_12.x;
        row[3 * self.index_2 + 1] = self.normal_12.y;
        row[3 * self.index_2 + 2] = (direct_orthogonal * self.point_2).dot(&self.normal_12);

        row
    }

    pub fn compute_jacobian_derivative(&self, states: &Vec<RigidBodyState>) -> RowDVector<Real> {
        let mut row = RowDVector::from(vec![0.0; 3 * states.len()]);
        row[3 * self.index_1 + 2] = states[self.index_1].velocity.z * self.point_1.dot(&self.normal_12);
        row[3 * self.index_2 + 2] = -states[self.index_2].velocity.z * self.point_2.dot(&self.normal_12);

        row
    }

    pub fn display(&self, rigid_bodies: &Vec<RigidBody>, context: Context, graphics: &mut impl Graphics, width: Scalar, height: Scalar) {
        let rigid_body = &rigid_bodies[self.index_1];
        let color = match self.kind {
            CollisionType::Penetration => [1.0, 0.0, 1.0, 0.5],
            CollisionType::Separation => [0.0, 1.0, 0.0, 0.5],
            CollisionType::Contact => [0.0, 0.0, 1.0, 0.5],
        };
        let point = Vector2::new(
            width as Real / 2.0 + (rigid_body.position.x + self.point_1.x) * 100.0,
            height as Real - (rigid_body.position.y + self.point_1.y) * 100.0
        );
        rectangle(
            color,
            [(point.x - 5.0) as Scalar, (point.y - 5.0) as Scalar, 10.0, 10.0],
            context.transform,
            graphics
        );

        let rigid_body = &rigid_bodies[self.index_2];
        let color = match self.kind {
            CollisionType::Penetration => [1.0, 0.0, 1.0, 0.5],
            CollisionType::Separation => [0.0, 1.0, 0.0, 0.5],
            CollisionType::Contact => [0.0, 0.0, 1.0, 0.5],
        };
        let point = Vector2::new(
            width as Real / 2.0 + (rigid_body.position.x + self.point_2.x) * 100.0,
            height as Real - (rigid_body.position.y + self.point_2.y) * 100.0
        );
        rectangle(
            color,
            [(point.x - 5.0) as Scalar, (point.y - 5.0) as Scalar, 10.0, 10.0],
            context.transform,
            graphics
        );
    }
}

pub fn compute_contact(index_1: Index, index_2: Index, rigid_bodies: &Vec<RigidBody>) -> Vec<Collision> {
    let rigid_body_1 = &rigid_bodies[index_1];
    let rigid_body_2 = &rigid_bodies[index_2];
    let mut manifolds: Vec<ContactManifold<(), ()>> = vec![];
    let transform = rigid_body_1.transform().inv_mul(rigid_body_2.transform());
    DefaultQueryDispatcher.contact_manifolds(
        &transform,
        rigid_body_1.shape(),
        rigid_body_2.shape(),
        0.0,
        &mut manifolds,
        &mut None
    ).expect("Collision non gérée !");

    let Some(contact) = manifolds.get(0) else {
        return vec![]
    };
    let result = find_deepest_contacts(contact);

    result.iter().map(|tracked_contact| {
        let direct_orthogonal = Rotation2::new(FRAC_PI_2);
        let rotation_1 = rigid_body_1.transform().rotation;
        let rotation_2 = rigid_body_2.transform().rotation;
        let point_1: Vector2<Real> = rotation_1 * tracked_contact.local_p1.coords;
        let point_2: Vector2<Real> = rotation_2 * tracked_contact.local_p2.coords;
        let normal_12: Vector2<Real> = rotation_1 * contact.local_n1;
        let velocity_1: Vector2<Real> = rigid_body_1.velocity.xy() + rigid_body_1.velocity.z * (direct_orthogonal * point_1);
        let velocity_2: Vector2<Real> = rigid_body_2.velocity.xy() + rigid_body_2.velocity.z * (direct_orthogonal * point_2);

        let relative_velocity_21 = (velocity_2 - velocity_1).dot(&normal_12);
        let kind;
        if relative_velocity_21.abs() < CONTACT_VELOCITY_EPS {
            kind = CollisionType::Contact;
        }
        else if relative_velocity_21 < 0.0 {
            kind = CollisionType::Penetration;
        }
        else {
            kind = CollisionType::Separation;
        }

        Collision {
            kind,
            index_1,
            index_2,
            point_1,
            point_2,
            normal_12
        }
    }).collect()
}

fn find_deepest_contacts(contact: &ContactManifold<(), ()>) -> Vec<&TrackedContact<()>> {
    let Some(mut deepest) = contact.points.get(0) else {
        return vec![];
    };
    let mut result = vec![deepest];

    for i in 1..contact.points.len() {
        if (contact.points[i].dist - deepest.dist).abs() < COLLISION_EPS {
            result.push(&contact.points[i]);
        }
        else if contact.points[i].dist < deepest.dist {
            result = vec![&contact.points[i]];
            deepest = &contact.points[i];
        }
    }

    result
}
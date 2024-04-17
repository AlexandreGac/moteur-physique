use std::f64::consts::FRAC_PI_2;
use nalgebra::{RowDVector, Vector2, Rotation2};
use parry2d_f64::math::Real;
use parry2d_f64::query::{ContactManifold, DefaultQueryDispatcher, PersistentQueryDispatcher, TrackedContact};
use parry2d_f64::utils::IsometryOpt;
use crate::rigid_body::RigidBody;
use crate::solver::Index;
use crate::utils::COLLISION_EPS;

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
        if relative_velocity_21 < 0.0 {
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

    result.iter().map(|tracked_contact| {let direct_orthogonal = Rotation2::new(FRAC_PI_2);
        let rotation_1 = contact.subshape_pos1.prepend_to(rigid_body_1.transform());
        let rotation_2 = contact.subshape_pos2.prepend_to(rigid_body_2.transform());
        let point_1: Vector2<Real> = rotation_1 * tracked_contact.local_p1.coords;
        let point_2: Vector2<Real> = rotation_2 * tracked_contact.local_p2.coords;
        let normal_12: Vector2<Real> = rotation_1 * contact.local_n1;
        let velocity_1: Vector2<Real> = rigid_body_1.velocity.xy() + rigid_body_1.velocity.z * (direct_orthogonal * point_1);
        let velocity_2: Vector2<Real> = rigid_body_2.velocity.xy() + rigid_body_2.velocity.z * (direct_orthogonal * point_2);

        let relative_velocity_21 = (velocity_2 - velocity_1).dot(&normal_12);
        let kind;
        if relative_velocity_21 < 0.0 {
            println!("Velocity : {relative_velocity_21}");
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
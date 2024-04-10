use nalgebra::{RowDVector, Point2, Vector2};
use parry2d::math::Real;
use parry2d::query::{ContactManifold, DefaultQueryDispatcher, PersistentQueryDispatcher, TrackedContact};
use crate::rigid_body::RigidBody;

#[derive(PartialEq, Clone, Debug)]
pub enum CollisionType {
    Penetration, Separation, Contact
}

#[derive(Clone, Debug)]
pub struct Collision {
    pub kind: CollisionType,
    pub local_point_1: Point2<Real>,
    pub local_point_2: Point2<Real>,
    pub local_normal_1: Vector2<Real>,
    pub local_normal_2: Vector2<Real>
}

impl Collision {

    pub fn update_collision_type(&mut self) {

    }

    pub fn get_restitution_coefficient(&self) -> Real {
        0.5
    }

    pub fn compute_jacobian(&self) -> RowDVector<Real> {
        RowDVector::<Real>::from(vec![1.0, 2.0, 3.0, 4.0, 5.0])
    }
}

pub fn compute_contact(rigid_body_1: &RigidBody, rigid_body_2: &RigidBody) -> Vec<Collision> {
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
        Collision {
            kind: CollisionType::Penetration,
            local_point_1: tracked_contact.local_p1,
            local_point_2: tracked_contact.local_p2,
            local_normal_1: contact.local_n1,
            local_normal_2: contact.local_n2
        }
    }).collect()
}

fn find_deepest_contacts(contact: &ContactManifold<(), ()>) -> Vec<&TrackedContact<()>> {
    let mut deepest = contact.points.get(0).expect("Contact non valide !");
    let mut result = vec![deepest];

    for i in 1..contact.points.len() {
        if contact.points[i].dist == deepest.dist {
            result.push(&contact.points[i]);
        }
        else if contact.points[i].dist < deepest.dist {
            result = vec![&contact.points[i]];
            deepest = &contact.points[i];
        }
    }

    result
}
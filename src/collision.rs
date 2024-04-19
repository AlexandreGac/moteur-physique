use std::f64::consts::FRAC_PI_2;
use nalgebra::{RowDVector, Vector2, Rotation2};
use parry2d_f64::math::Real;
use parry2d_f64::query::{ContactManifold, DefaultQueryDispatcher, PersistentQueryDispatcher, TrackedContact};
use parry2d_f64::shape::{FeatureId, ShapeType};
use piston_window::math::Scalar;
use piston_window::{Context, Graphics, rectangle};
use crate::rigid_body::{RigidBody, RigidBodyState};
use crate::solver::Index;
use crate::utils::{COLLISION_EPS, COLLISION_PREDICTION, CONTACT_VELOCITY_EPS, KINETIC_FRICTION, RESTITUTION_COEFFICIENT, STATIC_FRICTION};

#[derive(PartialEq, Clone, Debug)]
pub enum CollisionType {
    Penetration, Separation, Contact
}

#[derive(PartialEq, Clone, Debug)]
pub enum ContactType {
    Vertex, Edge, Curve
}

#[derive(Clone, Debug)]
pub struct Collision {
    pub kind: CollisionType,
    pub index_1: Index,
    pub index_2: Index,
    pub point_1: Vector2<Real>,
    pub point_2: Vector2<Real>,
    pub type_1: ContactType,
    pub type_2: ContactType,
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

    pub fn get_reaction_restitution_coefficient(&self, rigid_bodies: &Vec<RigidBody>, friction: bool, lagrangian: Option<&Real>, dt: Real) -> Real {
        if !friction {
            0.0
        }
        else {
            let lagrangian = *lagrangian.unwrap_or(&0.0);
            let direct_orthogonal = Rotation2::new(FRAC_PI_2);
            let rigid_body_1 = &rigid_bodies[self.index_1];
            let rigid_body_2 = &rigid_bodies[self.index_2];

            let orthogonal_point_1: Vector2<Real> = direct_orthogonal * self.point_1;
            let orthogonal_point_2: Vector2<Real> = direct_orthogonal * self.point_2;
            let velocity_1: Vector2<Real> = rigid_body_1.velocity.xy() + rigid_body_1.velocity.z * orthogonal_point_1;
            let velocity_2: Vector2<Real> = rigid_body_2.velocity.xy() + rigid_body_2.velocity.z * orthogonal_point_2;
            let relative_velocity_21: Vector2<Real> = velocity_2 - velocity_1;

            let tangent: Vector2<Real> = direct_orthogonal * self.normal_12;
            let tangential_velocity = relative_velocity_21.dot(&tangent).abs();

            let mass_term = rigid_body_1.inv_mass.x + rigid_body_2.inv_mass.x;
            let rotation_term_1 = orthogonal_point_1.dot(&tangent).powi(2) * rigid_body_1.inv_mass.z;
            let rotation_term_2 = orthogonal_point_2.dot(&tangent).powi(2) * rigid_body_2.inv_mass.z;
            let final_mass_term = mass_term + rotation_term_1 + rotation_term_2;

            if tangential_velocity < STATIC_FRICTION * lagrangian * dt * final_mass_term {
                return 0.0;
            }

            (tangential_velocity - KINETIC_FRICTION * lagrangian * dt * final_mass_term) / tangential_velocity
        }
    }

    pub fn get_impulse_restitution_coefficient(&self, rigid_bodies: &Vec<RigidBody>, friction: bool, lagrangian: Option<&Real>) -> Real {
        if !friction {
            RESTITUTION_COEFFICIENT
        }
        else {
            let lagrangian = *lagrangian.unwrap_or(&0.0);
            let direct_orthogonal = Rotation2::new(FRAC_PI_2);
            let rigid_body_1 = &rigid_bodies[self.index_1];
            let rigid_body_2 = &rigid_bodies[self.index_2];

            let orthogonal_point_1: Vector2<Real> = direct_orthogonal * self.point_1;
            let orthogonal_point_2: Vector2<Real> = direct_orthogonal * self.point_2;
            let velocity_1: Vector2<Real> = rigid_body_1.velocity.xy() + rigid_body_1.velocity.z * orthogonal_point_1;
            let velocity_2: Vector2<Real> = rigid_body_2.velocity.xy() + rigid_body_2.velocity.z * orthogonal_point_2;
            let relative_velocity_21: Vector2<Real> = velocity_2 - velocity_1;

            let tangent: Vector2<Real> = direct_orthogonal * self.normal_12;
            let tangential_velocity = relative_velocity_21.dot(&tangent).abs();

            let mass_term = rigid_body_1.inv_mass.x + rigid_body_2.inv_mass.x;
            let rotation_term_1 = orthogonal_point_1.dot(&tangent).powi(2) * rigid_body_1.inv_mass.z;
            let rotation_term_2 = orthogonal_point_2.dot(&tangent).powi(2) * rigid_body_2.inv_mass.z;
            let final_mass_term = mass_term + rotation_term_1 + rotation_term_2;

            if tangential_velocity < STATIC_FRICTION * lagrangian * final_mass_term {
                return 0.0;
            }

            -(tangential_velocity - KINETIC_FRICTION * lagrangian * final_mass_term) / tangential_velocity
        }
    }

    pub fn compute_jacobian(&self, rigid_bodies: &Vec<RigidBody>, friction: bool) -> RowDVector<Real> {
        let mut row = RowDVector::from(vec![0.0; 3 * rigid_bodies.len()]);
        let position_1 = rigid_bodies[self.index_1].position.xy();
        let position_2 = rigid_bodies[self.index_2].position.xy();
        let direct_orthogonal = Rotation2::new(FRAC_PI_2);
        let tangent = direct_orthogonal * self.normal_12;
        if !friction {
            row[3 * self.index_1] = -self.normal_12.x;
            row[3 * self.index_1 + 1] = -self.normal_12.y;
            match self.type_1 {
                ContactType::Edge => row[3 * self.index_1 + 2] = (position_2 + self.point_2 - position_1).dot(&tangent),
                _ => row[3 * self.index_1 + 2] = -(direct_orthogonal * self.point_1).dot(&self.normal_12)
            }

            row[3 * self.index_2] = self.normal_12.x;
            row[3 * self.index_2 + 1] = self.normal_12.y;
            match self.type_2 {
                ContactType::Edge => row[3 * self.index_2 + 2] = -(position_1 + self.point_1 - position_2).dot(&tangent),
                _ => row[3 * self.index_2 + 2] = (direct_orthogonal * self.point_2).dot(&self.normal_12)
            }
        }
        else {
            row[3 * self.index_1] = -tangent.x;
            row[3 * self.index_1 + 1] = -tangent.y;
            match self.type_1 {
                ContactType::Edge => row[3 * self.index_1 + 2] = -(position_2 + self.point_2 - position_1).dot(&self.normal_12),
                _ => row[3 * self.index_1 + 2] = -(direct_orthogonal * self.point_1).dot(&tangent)
            }

            row[3 * self.index_2] = tangent.x;
            row[3 * self.index_2 + 1] = tangent.y;
            match self.type_2 {
                ContactType::Edge => row[3 * self.index_2 + 2] = (position_1 + self.point_1 - position_2).dot(&self.normal_12),
                _ => row[3 * self.index_2 + 2] = (direct_orthogonal * self.point_2).dot(&tangent),
            }
        }

        row
    }

    pub fn compute_jacobian_derivative(&self, states: &Vec<RigidBodyState>) -> RowDVector<Real> {
        let mut row = RowDVector::from(vec![0.0; 3 * states.len()]);
        let position_1 = states[self.index_1].position.xy();
        let position_2 = states[self.index_2].position.xy();
        let velocity_1 = states[self.index_1].velocity.xy();
        let velocity_2 = states[self.index_2].velocity.xy();
        let direct_orthogonal = Rotation2::new(FRAC_PI_2);
        let tangent = direct_orthogonal * self.normal_12;
        match self.type_1 {
            ContactType::Vertex => row[3 * self.index_1 + 2] = states[self.index_1].velocity.z * self.point_1.dot(&self.normal_12),
            ContactType::Edge => row[3 * self.index_1 + 2] = 2.0 * (velocity_2 + states[self.index_2].velocity.z * (direct_orthogonal * self.point_2) - velocity_1).dot(&tangent) - states[self.index_1].velocity.z * (position_2 + self.point_2 - position_1).dot(&self.normal_12),
            ContactType::Curve => row[3 * self.index_1 + 2] = 0.0
        }

        match self.type_2 {
            ContactType::Vertex => row[3 * self.index_2 + 2] = -states[self.index_2].velocity.z * self.point_2.dot(&self.normal_12),
            ContactType::Edge => row[3 * self.index_2 + 2] = -2.0 * (velocity_1 + states[self.index_1].velocity.z * (direct_orthogonal * self.point_1) - velocity_2).dot(&tangent) + states[self.index_2].velocity.z * (position_1 + self.point_1 - position_2).dot(&self.normal_12),
            ContactType::Curve => row[3 * self.index_2 + 2] = 0.0
        }

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
            CollisionType::Contact  => [0.0, 0.0, 1.0, 0.5],
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
        COLLISION_PREDICTION,
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
        let kind = if relative_velocity_21.abs() < CONTACT_VELOCITY_EPS {
            CollisionType::Contact
        }
        else if relative_velocity_21 < 0.0 {
            CollisionType::Penetration
        }
        else {
            CollisionType::Separation
        };

        let type_1 = match tracked_contact.fid1.unpack() {
            FeatureId::Vertex(_) => ContactType::Vertex,
            FeatureId::Face(_) => if rigid_body_1.shape().shape_type() == ShapeType::Ball {
                ContactType::Curve
            } else { ContactType::Edge }
            _ => panic!("Type inconnu de collision !")
        };

        let type_2 = match tracked_contact.fid2.unpack() {
            FeatureId::Vertex(_) => ContactType::Vertex,
            FeatureId::Face(_) => if rigid_body_2.shape().shape_type() == ShapeType::Ball {
                ContactType::Curve
            } else { ContactType::Edge }
            _ => panic!("Type inconnu de collision !")
        };

        Collision {
            kind,
            index_1,
            index_2,
            point_1,
            point_2,
            type_1,
            type_2,
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

pub fn compute_friction_contacts(contacts: Vec<Collision>, lagrangians: &Vec<Real>) -> (Vec<Collision>, Vec<Real>) {
    if contacts.len() == 0 { return (vec![], vec![]) }
    let mut friction_contacts = vec![contacts[0].clone()];
    let mut friction_lagrangians = vec![lagrangians[0]];
    for i in 1..contacts.len() {
        if contacts[i - 1].index_1 == contacts[i].index_1 && contacts[i - 1].index_2 == contacts[i - 1].index_2 && contacts[i - 1].kind == contacts[i - 1].kind {
            friction_contacts.pop();
            let friction_contact = Collision {
                kind: contacts[i].kind.clone(),
                index_1: contacts[i].index_1,
                index_2: contacts[i].index_2,
                point_1: (contacts[i - 1].point_1 + contacts[i].point_1) / 2.0,
                point_2: (contacts[i - 1].point_2 + contacts[i].point_2) / 2.0,
                type_1: contacts[i].type_1.clone(),
                type_2: contacts[i].type_2.clone(),
                normal_12: contacts[i].normal_12
            };
            friction_contacts.push(friction_contact);

            friction_lagrangians.pop();
            friction_lagrangians.push(lagrangians[i - 1] + lagrangians[i]);
        }
        else {
            friction_contacts.push(contacts[i].clone());
            friction_lagrangians.push(lagrangians[i]);
        }
    }

    (friction_contacts, friction_lagrangians)
}
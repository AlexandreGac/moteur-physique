use crate::collision;
use crate::collision::{Collision, CollisionType};
use crate::rigid_body::{RigidBody, RigidBodyState};
use nalgebra::{DMatrix, DVector, Vector3};
use parry2d_f64::math::Real;
use piston_window::{Context, Graphics};
use piston_window::math::Scalar;
use crate::constraints::Constraint;
use crate::utils::GRAVITATIONAL_ACCELERATION;

pub struct Solver {
    dt: Real,
    world: Vec<RigidBody>,
    constraints: Vec<Constraint>,
    collision_infos: Vec<Collision>,
    contact_lagrangians: Vec<Real>,
    penetration_lagrangians: Vec<Real>,
    inverse_mass: DMatrix<Real>
}

impl Solver {

    pub fn simulate(&mut self) -> bool {
        let mut collisions = self.find_collisions();
        self.collision_infos = collisions.clone();
        let mut stop = false;
        loop {
            let contacts = collisions.iter().cloned()
                .filter(|collision| collision.kind == CollisionType::Contact)
                .collect::<Vec<_>>();

            let mut penetrations = collisions.iter().cloned()
                .filter(|collision| collision.kind == CollisionType::Penetration)
                .collect::<Vec<_>>();

            if !contacts.is_empty() {
                let impulses = self.compute_impulses(&contacts, &mut vec![], false);
                for i in 0..self.world.len() {
                    self.world[i].apply_impulse(impulses[i]);
                }
            }

            let impulses = self.compute_impulses(&vec![], &mut penetrations, false);
            for i in 0..impulses.len() {
                self.world[i].apply_impulse(impulses[i]);
            }

            if penetrations.is_empty() { break }

            let (mut friction_penetrations, friction_lagrangians) = collision::compute_friction_contacts(penetrations, &self.penetration_lagrangians);
            self.penetration_lagrangians = friction_lagrangians;
            let friction_impulses = self.compute_impulses(&vec![], &mut friction_penetrations, true);
            for i in 0..friction_impulses.len() {
                self.world[i].apply_impulse(friction_impulses[i]);
            }

            for collision in &mut collisions {
                collision.update_collision_type(&self.world);
            }

            stop = true;
        }

        let mut contacts = collisions.into_iter()
            .filter(|collision| collision.kind == CollisionType::Contact)
            .collect::<Vec<_>>();
        self.runge_kutta_4(&mut contacts);

        if !contacts.is_empty() {
            let (friction_contacts, friction_lagrangians) = collision::compute_friction_contacts(contacts, &self.contact_lagrangians);
            self.contact_lagrangians = friction_lagrangians;
            let friction_impulses = self.compute_impulses(&friction_contacts, &mut vec![], true);
            for i in 0..friction_impulses.len() {
                self.world[i].apply_impulse(friction_impulses[i]);
            }
        }

        println!("Énergie système : {}", self.compute_energy());
        stop
    }

    fn runge_kutta_4(&mut self, contacts: &mut Vec<Collision>) {
        loop {
            let states_1: Vec<RigidBodyState> = self.world.iter()
                .map(|x| x.get_state())
                .collect();
            let Some((forces_1, lagrangians_1)) = self.compute_forces(&states_1, contacts) else { continue };

            let states_2: Vec<RigidBodyState> = states_1.iter().enumerate()
                .map(|(i, x)| x.apply_physics_step(
                    states_1[i].velocity,
                    forces_1[i].component_mul(&self.world[i].inv_mass),
                    self.dt / 2.0
                )).collect();
            let Some((forces_2, lagrangians_2)) = self.compute_forces(&states_2, contacts) else { continue };

            let states_3: Vec<RigidBodyState> = states_1.iter().enumerate()
                .map(|(i, x)| x.apply_physics_step(
                    states_2[i].velocity,
                    forces_2[i].component_mul(&self.world[i].inv_mass),
                    self.dt / 2.0
                )).collect();
            let Some((forces_3, lagrangians_3)) = self.compute_forces(&states_3, contacts) else { continue };

            let states_4: Vec<RigidBodyState> = states_1.iter().enumerate()
                .map(|(i, x)| x.apply_physics_step(
                    states_3[i].velocity,
                    forces_3[i].component_mul(&self.world[i].inv_mass),
                    self.dt
                )).collect();
            let Some((forces_4, lagrangians_4)) = self.compute_forces(&states_4, contacts) else { continue };

            for i in 0..self.world.len() {
                let velocity = states_1[i].velocity + 2.0 * (states_2[i].velocity + states_3[i].velocity) + states_4[i].velocity;
                let acceleration = self.world[i].inv_mass.component_mul(&(forces_1[i] + 2.0 * (forces_2[i] + forces_3[i]) + forces_4[i]));
                self.world[i].apply_physics(velocity, acceleration, self.dt / 6.0);
            }

            self.contact_lagrangians = vec![];
            for i in 0..contacts.len() {
                self.contact_lagrangians.push((lagrangians_1[i] + 2.0 * (lagrangians_2[i] + lagrangians_3[i]) + lagrangians_4[i]) / 6.0);
            }

            break
        }
    }

    fn find_collisions(&self) -> Vec<Collision> {
        let n = self.world.len();
        let mut contacts = vec![];
        for i in 0..n {
            for j in (i + 1)..n {
                contacts.push(collision::compute_contact(i, j, &self.world));
            }
        }

        contacts.into_iter().flatten().collect()
    }

    fn compute_impulses(&mut self, contacts: &Vec<Collision>, penetrations: &mut Vec<Collision>, friction: bool) -> Vec<Vector3<Real>> {
        let states: Vec<RigidBodyState> = self.world.iter().map(|x| x.get_state()).collect();
        loop {
            if contacts.len() == 0 && penetrations.len() == 0 { return vec![] }

            let velocities: DVector<Real> = DVector::<Real>::from(
                self.world.iter()
                    .map(|x| x.velocity.iter().cloned())
                    .flatten()
                    .collect::<Vec<_>>()
            );
            let jacobian: DMatrix<Real> = DMatrix::<Real>::from_rows(
                vec![
                    self.constraints.iter()
                        .map(|x| x.compute_jacobian(&states))
                        .collect::<Vec<_>>(),
                    contacts.iter()
                        .map(|x| x.compute_jacobian(&self.world, friction))
                        .collect::<Vec<_>>(),
                    penetrations.iter()
                        .map(|x| x.compute_jacobian(&self.world, friction))
                        .collect::<Vec<_>>()
                ]
                    .into_iter().flatten().collect::<Vec<_>>().as_slice()
            );
            let bias: DMatrix<Real> = DMatrix::<Real>::from_diagonal(
                &vec![
                    vec![1.0; self.constraints.len()],
                    contacts.iter().enumerate()
                        .map(|(i, x)| 1.0 - x.get_reaction_restitution_coefficient(&self.world, friction, self.contact_lagrangians.get(i), self.dt))
                        .collect::<Vec<_>>(),
                    penetrations.iter().enumerate()
                        .map(|(i, x)| 1.0 + x.get_impulse_restitution_coefficient(&self.world, friction, self.penetration_lagrangians.get(i)))
                        .collect::<Vec<_>>()
                ]
                    .into_iter().flatten().collect::<Vec<_>>().into()
            );

            let a = jacobian.clone() * &self.inverse_mass * jacobian.transpose();
            let b = -bias * jacobian.clone() * velocities;
            let lagrangian = a.svd(true, true).solve(&b, 1e-15).expect("La résolution a échoué !");

            if !friction {
                let mut restart = false;
                for i in (0..penetrations.len()).rev() {
                    if lagrangian[self.constraints.len() + contacts.len() + i] < 0.0 {
                        penetrations.remove(i);
                        restart = true;
                    }
                }
                if restart { continue }

                self.penetration_lagrangians = vec![];
                for i in 0..penetrations.len() {
                    self.penetration_lagrangians.push(lagrangian[self.constraints.len() + contacts.len() + i]);
                }
            }

            let constraint_impulses = jacobian.transpose() * lagrangian;

            let length = constraint_impulses.len() / 3;
            let mut impulses = vec![];
            for i in 0..length {
                impulses.push(Vector3::new(
                    constraint_impulses[3 * i],
                    constraint_impulses[3 * i + 1],
                    constraint_impulses[3 * i + 2]
                ));
            }

            return impulses;
        }
    }

    fn compute_forces(&mut self, states: &Vec<RigidBodyState>, contacts: &mut Vec<Collision>) -> Option<(Vec<Vector3<Real>>, Vec<Real>)> {
        let mut forces = vec![Vector3::<Real>::new(0.0, 0.0, 0.0); states.len()];
        for i in 0..states.len() {
            let inv_mass = self.world[i].inv_mass.x;
            if inv_mass > 0.0 {
                forces[i] = Vector3::<Real>::new(0.0, -GRAVITATIONAL_ACCELERATION / inv_mass, 0.0);
            }
        }

        if self.constraints.len() == 0 && contacts.len() == 0 { return Some((forces, vec![])) }

        let velocities: DVector<Real> = DVector::<Real>::from(
            states.iter()
                .map(|x| x.velocity.iter().cloned())
                .flatten()
                .collect::<Vec<_>>()
        );
        let external_forces: DVector<Real> = DVector::<Real>::from(
            forces.iter()
                .map(|x| x.iter().cloned())
                .flatten()
                .collect::<Vec<_>>()
        );
        let jacobian: DMatrix<Real> = DMatrix::<Real>::from_rows(
            vec![
                self.constraints.iter()
                    .map(|x| x.compute_jacobian(&states))
                    .collect::<Vec<_>>(),
                contacts.iter()
                    .map(|x| x.compute_jacobian(&self.world, false))
                    .collect::<Vec<_>>()
            ]
                .into_iter().flatten().collect::<Vec<_>>().as_slice()
        );
        let jacobian_derivative: DMatrix<Real> = DMatrix::<Real>::from_rows(
            vec![
                self.constraints.iter()
                    .map(|x| x.compute_jacobian_derivative(&states))
                    .collect::<Vec<_>>(),
                contacts.iter()
                    .map(|x| x.compute_jacobian_derivative(&states))
                    .collect::<Vec<_>>()
            ]
                .into_iter().flatten().collect::<Vec<_>>().as_slice()
        );

        let a = jacobian.clone() * &self.inverse_mass * jacobian.transpose();
        let b = -jacobian_derivative * velocities
            - jacobian.clone() * &self.inverse_mass * external_forces;
        let lagrangian = a.svd(true, true).solve(&b, 1e-15).expect("La résolution a échoué !");

        let mut restart = false;
        for i in (0..contacts.len()).rev() {
            if lagrangian[self.constraints.len() + i] < 0.0 {
                contacts.remove(i);
                restart = true;
            }
        }
        if restart { return None; }

        let mut lagrangians = vec![];
        for i in 0..contacts.len() {
            lagrangians.push(lagrangian[self.constraints.len() + i]);
        }

        let constraint_forces = jacobian.transpose() * lagrangian;

        let length = states.len();
        let mut result_forces = vec![];
        for i in 0..length {
            result_forces.push(Vector3::new(
                constraint_forces[3 * i] + forces[i].x,
                constraint_forces[3 * i + 1] + forces[i].y,
                constraint_forces[3 * i + 2] + forces[i].z
            ));
        }

        Some((result_forces, lagrangians))
    }

    pub fn compute_energy(&self) -> Real {
        self.world.iter().filter(|x| x.inv_mass.x > 0.0).map(|x| x.get_energy()).sum()
    }

    pub fn display(&self, context: Context, graphics: &mut impl Graphics, width: Scalar, height: Scalar) {
        self.world.iter().for_each(|x| x.display(context, graphics, width, height));
        self.constraints.iter().for_each(|x| x.display(&self.world, context, graphics, width, height));
        self.collision_infos.iter().for_each(|x| x.display(&self.world, context, graphics, width, height));
    }
}

pub type Index = usize;

pub struct SolverBuilder {
    world: Vec<RigidBody>,
    constraints: Vec<Constraint>
}

impl SolverBuilder {

    pub fn new() -> Self {
        Self {
            world: vec![],
            constraints: vec![]
        }
    }

    pub fn add_rigid_body(mut self, rigid_body: RigidBody) -> Self {
        self.world.push(rigid_body);
        self
    }

    pub fn add_constraint(mut self, constraint: Constraint) -> Self {
        self.constraints.push(constraint);
        self
    }

    pub fn build(self) -> Solver {
        let dt = 1.0 / 600.0;
        let inverse_mass: DMatrix<Real> = DMatrix::<Real>::from_diagonal(
            &self.world.iter()
                .map(|x| x.inv_mass.iter())
                .flatten()
                .cloned()
                .collect::<Vec<_>>()
                .into()
        );

        Solver {
            dt,
            world: self.world,
            constraints: self.constraints,
            collision_infos: vec![],
            contact_lagrangians: vec![],
            penetration_lagrangians: vec![],
            inverse_mass
        }
    }
}
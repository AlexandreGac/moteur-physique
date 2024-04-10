use crate::collision;
use crate::collision::{Collision, CollisionType};
use crate::rigid_body::{RigidBody, RigidBodyState};
use nalgebra::{DMatrix, DVector, Vector3};
use parry2d::math::Real;
use crate::constraints::Constraint;

pub struct Solver {
    dt: Real,
    world: Vec<RigidBody>,
    constraints: Vec<Constraint>,
    inverse_mass_properties: DMatrix<Real>
}

impl Solver {

    pub fn simulate(&mut self) {
        // Replace with a loop to prevent code repetition ?
        let mut collisions = self.find_collisions();
        let mut penetrations = collisions.iter().cloned()
            .filter(|collision| collision.kind == CollisionType::Penetration)
            .collect::<Vec<_>>();

        while !penetrations.is_empty() {
            let impulses = self.compute_impulses(vec![], penetrations);
            for i in 0..self.world.len() {
                self.world[i].apply_impulse(impulses[i]);
            }
            for collision in &mut collisions {
                collision.update_collision_type();
            }
            penetrations = collisions.iter().cloned()
                .filter(|collision| collision.kind == CollisionType::Penetration)
                .collect::<Vec<_>>();
        }

        self.runge_kutta_4();
    }

    fn runge_kutta_4(&mut self) {
        let states_1: Vec<RigidBodyState> = self.world.iter()
            .map(|x| x.get_state())
            .collect();
        let forces_1 = self.compute_forces(&states_1, vec![]);

        let states_2: Vec<RigidBodyState> = states_1.iter().enumerate()
            .map(|(i, x)| x.apply_physics_step(
                states_1[i].velocity,
                forces_1[i].component_mul(&self.world[i].inv_mass),
                self.dt / 2.0
            )).collect();
        let forces_2 = self.compute_forces(&states_2, vec![]);

        let states_3: Vec<RigidBodyState> = states_1.iter().enumerate()
            .map(|(i, x)| x.apply_physics_step(
                states_2[i].velocity,
                forces_2[i].component_mul(&self.world[i].inv_mass),
                self.dt / 2.0
            )).collect();
        let forces_3 = self.compute_forces(&states_3, vec![]);

        let states_4: Vec<RigidBodyState> = states_1.iter().enumerate()
            .map(|(i, x)| x.apply_physics_step(
                states_3[i].velocity,
                forces_3[i].component_mul(&self.world[i].inv_mass),
                self.dt
            )).collect();
        let forces_4 = self.compute_forces(&states_4, vec![]);

        for i in 0..self.world.len() {
            let velocity = states_1[i].velocity + 2.0 * (states_2[i].velocity + states_3[i].velocity) + states_4[i].velocity;
            let acceleration = self.world[i].inv_mass.component_mul(&(forces_1[i] + 2.0 * (forces_2[i] + forces_3[i]) + forces_4[i]));
            self.world[i].apply_physics(velocity, acceleration, self.dt);
        }
    }

    fn find_collisions(&self) -> Vec<Collision> {
        let n = self.world.len();
        let mut contacts = vec![];
        for i in 0..n {
            for j in (i + 1)..n {
                contacts.push(collision::compute_contact(
                    &self.world[i],
                    &self.world[j])
                );
            }
        }

        contacts.into_iter().flatten().collect()
    }

    fn compute_impulses(&self, constraints: Vec<Constraint>, penetrations: Vec<Collision>) -> Vec<Vector3<Real>> {
        let velocities: DVector<Real> = DVector::<Real>::from(
            self.world.iter()
                .map(|x| x.velocity.iter().cloned())
                .flatten()
                .collect::<Vec<_>>()
        );
        let jacobian: DMatrix<Real> = DMatrix::<Real>::from_rows(
            vec![
                constraints.iter()
                    .map(|x| x.compute_jacobian())
                    .collect::<Vec<_>>(),
                penetrations.iter()
                    .map(|x| x.compute_jacobian())
                    .collect::<Vec<_>>()
            ]
                .into_iter().flatten().collect::<Vec<_>>().as_slice()
        );
        let bias: DMatrix<Real> = DMatrix::<Real>::from_diagonal(
            &penetrations.iter()
                .map(|x| 1.0 + x.get_restitution_coefficient())
                .collect::<Vec<_>>()
                .into()
        );

        let a = jacobian.clone() * &self.inverse_mass_properties * jacobian.transpose();
        let b = -bias * jacobian.clone() * velocities;
        let lagrangian = a.lu().solve(&b).expect("La résolution a échoué !");
        let constraint_impulses = jacobian.transpose() * lagrangian;

        let length = constraint_impulses.len() / 3;
        let mut impulses = vec![];
        for i in 0..length {
            impulses.push(Vector3::new(
                constraint_impulses[i],
                constraint_impulses[i + 1],
                constraint_impulses[i + 2]
            ));
        }

        impulses
    }

    fn compute_forces(&self, states: &Vec<RigidBodyState>, constraints: Vec<Constraint>) -> Vec<Vector3<Real>> {
        let velocities: DVector<Real> = DVector::<Real>::from(
            states.iter()
                .map(|x| x.velocity.iter().cloned())
                .flatten()
                .collect::<Vec<_>>()
        );
        let external_forces: DVector<Real> = DVector::<Real>::from(
            vec![0.0; states.len()]
        );
        let jacobian: DMatrix<Real> = DMatrix::<Real>::from_rows(
            constraints.iter()
                .map(|x| x.compute_jacobian())
                .collect::<Vec<_>>()
                .as_slice()
        );
        let jacobian_derivative: DMatrix<Real> = DMatrix::<Real>::from_rows(
            constraints.iter()
                .map(|x| x.compute_jacobian_derivative())
                .collect::<Vec<_>>()
                .as_slice()
        );

        let a = jacobian.clone() * &self.inverse_mass_properties * jacobian.transpose();
        let b = -jacobian_derivative * velocities
            - jacobian.clone() * &self.inverse_mass_properties * external_forces;
        let lagrangian = a.lu().solve(&b).expect("La résolution a échoué !");
        let constraint_forces = jacobian.transpose() * lagrangian;

        let length = constraint_forces.len() / 3;
        let mut forces = vec![];
        for i in 0..length {
            forces.push(Vector3::new(
                constraint_forces[i],
                constraint_forces[i + 1],
                constraint_forces[i + 2]
            ));
        }

        forces
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
        let dt = 1.0 / 60.0;
        let inverse_mass_properties: DMatrix<Real> = DMatrix::<Real>::from_diagonal(
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
            inverse_mass_properties
        }
    }
}
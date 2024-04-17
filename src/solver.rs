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
    inverse_mass: DMatrix<Real>
}

impl Solver {

    pub fn simulate(&mut self) {
        let mut collisions = self.find_collisions();
        loop {
            let penetrations = collisions.iter().cloned()
                .filter(|collision| collision.kind == CollisionType::Penetration)
                .collect::<Vec<_>>();

            if penetrations.is_empty() { break }

            let impulses = self.compute_impulses(penetrations);
            for i in 0..self.world.len() {
                self.world[i].apply_impulse(impulses[i]);
            }
            for collision in &mut collisions {
                collision.update_collision_type(&self.world);
            }
        }

        self.runge_kutta_4();
        println!("Énergie système : {}", self.compute_energy());
    }

    fn runge_kutta_4(&mut self) {
        let states_1: Vec<RigidBodyState> = self.world.iter()
            .map(|x| x.get_state())
            .collect();
        let forces_1 = self.compute_forces(&states_1);

        let states_2: Vec<RigidBodyState> = states_1.iter().enumerate()
            .map(|(i, x)| x.apply_physics_step(
                states_1[i].velocity,
                forces_1[i].component_mul(&self.world[i].inv_mass),
                self.dt / 2.0
            )).collect();
        let forces_2 = self.compute_forces(&states_2);

        let states_3: Vec<RigidBodyState> = states_1.iter().enumerate()
            .map(|(i, x)| x.apply_physics_step(
                states_2[i].velocity,
                forces_2[i].component_mul(&self.world[i].inv_mass),
                self.dt / 2.0
            )).collect();
        let forces_3 = self.compute_forces(&states_3);

        let states_4: Vec<RigidBodyState> = states_1.iter().enumerate()
            .map(|(i, x)| x.apply_physics_step(
                states_3[i].velocity,
                forces_3[i].component_mul(&self.world[i].inv_mass),
                self.dt
            )).collect();
        let forces_4 = self.compute_forces(&states_4);

        for i in 0..self.world.len() {
            let velocity = states_1[i].velocity + 2.0 * (states_2[i].velocity + states_3[i].velocity) + states_4[i].velocity;
            let acceleration = self.world[i].inv_mass.component_mul(&(forces_1[i] + 2.0 * (forces_2[i] + forces_3[i]) + forces_4[i]));
            self.world[i].apply_physics(velocity, acceleration, self.dt / 6.0);
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

    fn compute_impulses(&self, mut penetrations: Vec<Collision>) -> Vec<Vector3<Real>> {
        let states: Vec<RigidBodyState> = self.world.iter().map(|x| x.get_state()).collect();
        loop {
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
                    penetrations.iter()
                        .map(|x| x.compute_jacobian(&self.world))
                        .collect::<Vec<_>>()
                ]
                    .into_iter().flatten().collect::<Vec<_>>().as_slice()
            );
            let bias: DMatrix<Real> = DMatrix::<Real>::from_diagonal(
                &vec![
                    vec![1.0; self.constraints.len()],
                    penetrations.iter()
                        .map(|x| 1.0 + x.get_restitution_coefficient())
                        .collect::<Vec<_>>()
                ]
                    .into_iter().flatten().collect::<Vec<_>>().into()
            );

            let a = jacobian.clone() * &self.inverse_mass * jacobian.transpose();
            let b = -bias * jacobian.clone() * velocities;
            let lagrangian = a.svd(true, true).solve(&b, 1e-15).expect("La résolution a échoué !");

            let mut restart = false;
            for i in (0..penetrations.len()).rev() {
                if lagrangian[self.constraints.len() + i] < 0.0 {
                    penetrations.remove(i);
                    restart = true;
                }
            }
            if restart { continue }

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

    fn compute_forces(&self, states: &Vec<RigidBodyState>) -> Vec<Vector3<Real>> {
        let mut forces = vec![Vector3::<Real>::new(0.0, 0.0, 0.0); states.len()];
        for i in 0..states.len() {
            let inv_mass = self.world[i].inv_mass.x;
            if inv_mass > 0.0 {
                forces[i] = Vector3::<Real>::new(0.0, -GRAVITATIONAL_ACCELERATION / inv_mass, 0.0);
            }
        }

        if self.constraints.len() == 0 { return forces }

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
            self.constraints.iter()
                .map(|x| x.compute_jacobian(states))
                .collect::<Vec<_>>()
                .as_slice()
        );
        let jacobian_derivative: DMatrix<Real> = DMatrix::<Real>::from_rows(
            self.constraints.iter()
                .map(|x| x.compute_jacobian_derivative(states))
                .collect::<Vec<_>>()
                .as_slice()
        );

        let a = jacobian.clone() * &self.inverse_mass * jacobian.transpose();
        let b = -jacobian_derivative * velocities
            - jacobian.clone() * &self.inverse_mass * external_forces;
        let lagrangian = a.lu().solve(&b).expect("La résolution a échoué !");
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

        result_forces
    }

    pub fn compute_energy(&self) -> Real {
        self.world.iter().filter(|x| x.inv_mass.x > 0.0).map(|x| x.get_energy()).sum()
    }

    pub fn display(&self, context: Context, graphics: &mut impl Graphics, width: Scalar, height: Scalar) {
        self.world.iter().for_each(|x| x.display(context, graphics, width, height));
        self.constraints.iter().for_each(|x| x.display(&self.world, context, graphics, width, height));
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
            inverse_mass
        }
    }
}
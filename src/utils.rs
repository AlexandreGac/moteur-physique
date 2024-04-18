use parry2d_f64::math::Real;

pub const GRAVITATIONAL_ACCELERATION: Real = 9.81;
pub const COLLISION_EPS: Real = 1e-6;
pub const CONTACT_VELOCITY_EPS: Real = 2e-2;  // Valeur empirique, liée au dt
pub const RESTITUTION_COEFFICIENT: Real = 0.5;
pub const STATIC_FRICTION: Real = 0.31;
pub const KINETIC_FRICTION: Real = 0.23;